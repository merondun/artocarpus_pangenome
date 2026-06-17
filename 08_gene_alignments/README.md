# 08_gene_alignments/

Gene-based whole genome alignment using JCVI. 

Outputs:

![jcvi](/imgs/whole_genome_alignment_gene.png)



___



## Prep

* Extracts longest isoforms per gene, generates CDS/PEP FASTA files, validates start codons, and prepares JCVI‑compatible BED/CDS files for downstream gene‑alignment workflows.

First, extract cds and essential files for single longest transcript per gene:

```bash
#!/bin/bash

#SBATCH --time=1-00:00:00   
#SBATCH --cpus-per-task=2
#SBATCH --mem=8Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

WD=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/annotation
mkdir -p ${WD}/only_longest_transcript_per_gene
cd ${WD}

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <sample>"
    exit 1
fi
# module load miniconda
# source activate isoseq_ann

SAMPLE=$1

TARGET=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/${SAMPLE}.fa

agat_sp_keep_longest_isoform.pl --gff gtfs/${SAMPLE}.gtf -o only_longest_transcript_per_gene/${SAMPLE}.longest_transcript_per_gene.gtf
gffread only_longest_transcript_per_gene/${SAMPLE}.longest_transcript_per_gene.gtf -o only_longest_transcript_per_gene/${SAMPLE}.gff3 --keep-genes -O -g ${TARGET} -w only_longest_transcript_per_gene/${SAMPLE}.fa
gffread only_longest_transcript_per_gene/${SAMPLE}.gff3 \
    -g ${TARGET} \
    -x only_longest_transcript_per_gene/${SAMPLE}.cds.fa

gffread only_longest_transcript_per_gene/${SAMPLE}.gff3 \
    -g ${TARGET} \
    -y only_longest_transcript_per_gene/${SAMPLE}.pep.fa
```

And for camansi / mari:

```bash
for SAMPLE in HART063 HART067; do 
    TARGET=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/${SAMPLE}.fa
    gffread only_longest_transcript_per_gene/${SAMPLE}.gff3 \
        -g ${TARGET} \
        -x only_longest_transcript_per_gene/${SAMPLE}.cds.fa

    gffread only_longest_transcript_per_gene/${SAMPLE}.gff3 \
        -g ${TARGET} \
        -y only_longest_transcript_per_gene/${SAMPLE}.pep.fa
done 
```

Submit serial via sample:

```bash
cat Samples.list | xargs -I {} sbatch -J ann_{} 01_Extra_Single_Transcript.sh {} 
```

Confirm starts:

```bash
echo -e "Sample\tTotal\tATG_start\tM_start" > start_summary.tsv

for CDS in *.cds.fa; do
    SAMPLE=${CDS%.cds.fa}
    PEP=${SAMPLE}.pep.fa

    TOTAL=$(grep -c "^>" "$CDS")

    ATGstart=$(seqkit seq -w 0 "$CDS" | \
        awk '/^>/{getline seq; if (substr(seq,1,3)=="ATG") c++} END{print c+0}')

    Mstart=$(seqkit seq -w 0 "$PEP" | \
        awk '/^>/{getline seq; if (substr(seq,1,1)=="M") c++} END{print c+0}')

    echo -e "${SAMPLE}\t${TOTAL}\t${ATGstart}\t${Mstart}" >> start_summary.tsv
done
```

Output:

```bash
Sample  Total   ATG_start       M_start
H6      23463   23433   23433
HART001 23890   23858   23858
HART030 23790   23762   23762
HART032 23678   23645   23645
HART033 23669   23637   23637
HART038 23689   23661   23661
HART046 23516   23487   23487
HART049 23150   23123   23123
HART050 23818   23791   23791
HART053 23540   23506   23506
HART069 23652   23622   23622
ZZ3     23343   23314   23314
ZZ7     23261   23232   23232
ZZ9     23478   23446   23446
HART063 27770   27746   27746
HART067 26601   26052   26052
```

Copy over:

```bash
WD=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/annotation/only_longest_transcript_per_gene
cd ${WD}
mkdir -p ../jcvi
for i in $(cat ../Samples.list); do 
	echo "Working on ${i}"
	cp ${i}.cds.fa ../jcvi/${i}.cds.fa 	
	cp ${i}.longest_transcript_per_gene.gtf ../jcvi/${i}.gff3

done
```

Convert to bed with jcvi:

```bash
WD=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/annotation/jcvi
cd ${WD}
# source activate jcvi
for i in $(cat  ../Samples.list); do 
	echo "Working on ${i}"
	python -m jcvi.formats.gff bed --type=transcript --key=ID ${i}.gff3 -o ${i}.bed
	python -m jcvi.formats.fasta format ${i}.cds.fa ${i}.cds
done 
```



## Alignments

* Runs pairwise genome synteny with JCVI (orthologs, dotplots, filtered anchors) and annotates inter‑chromosomal anchors.

Synteny with JCVI:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=10
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate jcvi
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <first_pair> <second_pair>"
    exit 1
fi

reference=$1  #reference
target=$2 #target
echo "Working on $target aligned to $reference"

WD=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/annotation/jcvi
cd ${WD}

rm ${reference}*last* ${reference}*lifted* ${reference}*pdf

# Align and extract. Note the formatting must be EXACT with file names, cscore = RBH according to jcvi documentation
python -m jcvi.compara.catalog ortholog ${reference} ${target} --cscore=.99 
#--no_strip_names 

# for dotplots
magick -density 300 ${reference}.${target}.pdf ../output_jcvi/${reference}_${target}.png

# for filtering anchors 
python -m jcvi.compara.synteny screen --minsize=10 --simple ${reference}.${target}.anchors ${reference}.${target}.anchors.minsize10

```

Submit:

```bash
cat Pairs.list | xargs -L 1 sbatch 01_JCVI.sh
```

For inter-chromosomal colors:

```bash
#!/usr/bin/env bash

# Usage:
#   ./label_interchrom.sh Pairs.list
#
# inputs:
#   ${ref}.bed
#   ${tgt}.bed
#   ${ref}.${tgt}.anchors.simple
#
# Produces:
#   ${ref}.${tgt}.anchors.cols.simple   (annotated with g* for inter-chr)

pairs_file="$1"

while read -r ref tgt; do
    echo "Processing pair: $ref  $tgt"

    ref_bed="${ref}.bed"
    tgt_bed="${tgt}.bed"
    anchor_file="${ref}.${tgt}.anchors.simple"
    out_file="${ref}.${tgt}.anchors.cols.simple"

    if [[ ! -f "$ref_bed" || ! -f "$tgt_bed" || ! -f "$anchor_file" ]]; then
        echo "  Missing files for pair ($ref, $tgt). Skipping."
        continue
    fi

    # Build chromosome lookup tables
    awk '{chr[$4]=$1} END {for (k in chr) print k, chr[k]}' "$ref_bed" > ref.chr.map
    awk '{chr[$4]=$1} END {for (k in chr) print k, chr[k]}' "$tgt_bed" > tgt.chr.map

    # Annotate anchors
    awk '
        BEGIN {
            while (getline <"ref.chr.map" > 0) { refchr[$1]=$2 }
            while (getline <"tgt.chr.map" > 0) { tgtchr[$1]=$2 }
        }
        {
            r1=$1; r2=$2; t1=$3; t2=$4;
            inter = (refchr[r1] != tgtchr[t1] || refchr[r2] != tgtchr[t2])

            if (inter)
                print "g*"$0
            else
                print $0
        }
    ' "$anchor_file" > "$out_file"

    echo "  -> wrote $out_file"

done < "$pairs_file"

# Cleanup temporary maps
```



## Plot 

Plot with JCVI. 

Create a layout file:

```
#	y,	xstart,	xend,	rotation,	color,	label,	va,	bed
	0.94,	0.2,	0.95,	0,	,	HART001,	top,	HART001.bed
	0.88,	0.2,	0.95,	0,	,	HART030,	top,	HART030.bed
	0.82,	0.2,	0.95,	0,	,	HART050,	top,	HART050.bed
	0.76,	0.2,	0.95,	0,	,	HART069,	top,	HART069.bed
	0.71,	0.2,	0.95,	0,	,	HART032,	top,	HART032.bed
	0.65,	0.2,	0.95,	0,	,	HART033,	top,	HART033.bed
	0.59,	0.2,	0.95,	0,	,	HART038,	top,	HART038.bed
	0.53,	0.2,	0.95,	0,	,	ZZ3,	top,	ZZ3.bed
	0.47,	0.2,	0.95,	0,	,	ZZ7,	top,	ZZ7.bed
	0.41,	0.2,	0.95,	0,	,	ZZ9,	top,	ZZ9.bed
	0.35,	0.2,	0.95,	0,	,	HART046,	top,	HART046.bed
	0.29,	0.2,	0.95,	0,	,	HART049,	top,	HART049.bed
	0.24,	0.2,	0.95,	0,	,	HART053,	top,	HART053.bed
	0.18,	0.2,	0.95,	0,	,	H6,	top,	H6.bed
	0.12,	0.2,	0.95,	0,	,	HART063,	top,	HART063.bed
	0.06,	0.2,	0.95,	0,	,	HART067,	bottom,	HART067.bed
#	edges							
e,	0,	1,	HART001.HART030.anchors.cols.simple					
e,	1,	2,	HART030.HART050.anchors.cols.simple					
e,	2,	3,	HART050.HART069.anchors.cols.simple					
e,	3,	4,	HART069.HART032.anchors.cols.simple					
e,	4,	5,	HART032.HART033.anchors.cols.simple					
e,	5,	6,	HART033.HART038.anchors.cols.simple					
e,	6,	7,	HART038.ZZ3.anchors.cols.simple					
e,	7,	8,	ZZ3.ZZ7.anchors.cols.simple					
e,	8,	9,	ZZ7.ZZ9.anchors.cols.simple					
e,	9,	10,	ZZ9.HART046.anchors.cols.simple					
e,	10,	11,	HART046.HART049.anchors.cols.simple					
e,	11,	12,	HART049.HART053.anchors.cols.simple					
e,	12,	13,	HART053.H6.anchors.cols.simple					
e,	13,	14,	H6.HART063.anchors.cols.simple					
e,	14,	15,	HART063.HART067.anchors.cols.simple					
```

And the chrs.txt file. I took this ordering from the chromsyn ordering. 

```
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
```

And create karyoplot:

```bash
mkdir -p ../output_jcvi
python -m jcvi.graphics.karyotype \
	--format png --font Arial --seed 1 \
	-o ../output_jcvi/20260616_gene_alignments_interchr_min10.pdf \
	chrs.txt chr_layout.txt --basepair
```