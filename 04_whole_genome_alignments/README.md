# 04_whole_genome_alignments/

- confirm broad-scale chromosomal synteny and orientation across all accessions using whole‑genome dotplots (HART001 as reference) 

- correct a putative misassembly: invert Chr07 in HART001 based on consistent dotplot inversion + Hi‑C support, ensuring all genomes share a comparable Chr07 direction before pangenome alignment 



Outputs, chromosomal synteny among *Artocarpus altilis* : 

![dotplots](/imgs/dotplots.png)



BUSCO-based whole genome alignments : 

![chromsyn](/imgs/busco_alignments.png)

___



## Dotplots against HART001

Using HART001 as the reference, produce dotplots:

```bash
#!/bin/bash

#SBATCH --time=1-00:00:00    
#SBATCH --cpus-per-task=10
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

SAMPLE=$1

# module load minicionda
# source activate puzzler200
module load apptainer

WD=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/dotplots
GENOMES=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked
REF=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/primary_asm/init/HART001.Chr07sub.fa

mkdir -p pafs

echo "WGA for ${SAMPLE}"
mashmap -r ${REF} -q ${GENOMES}/${SAMPLE}.chr.fa -t 10 -s 10000 --perc_identity 95 -o pafs/${SAMPLE}.paf 2> pafs/${SAMPLE}.mashmap.log
Rscript ~/apptainer/paf2dotplot.R pafs/${SAMPLE}.paf -r 1e6 -m 1e4 -p 4 -c 1 --sort-by-refid --identity-lower-color 100
```

Submit: 

```bash
cat Samples.list  | xargs -I {} sbatch -J wga_{} 01_WGA_HART001Ref.sh {}
zip wga_dotplots_hart001.zip pafs/*pdf
```

## BUSCO-level Synteny

This workflow runs chromsyn to place BUSCO genes onto chromosomes and summarize synteny using BUSCO anchors. It generates plotting inputs (BUSCO tables, telomere tracks, and repeat/telomere-window scores), merges them into a chromsyn fig. 

Files:

```
cat Fastas.list 
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/H6.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/HART001.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/HART030.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/HART032.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/HART033.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/HART038.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/HART046.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/HART049.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/HART050.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/HART053.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/HART063.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/HART067.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/HART069.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/ZZ3.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/ZZ7.chr.fa
/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/ZZ9.chr.fa
```

Run the chromsyn workflow:

```bash
#!/bin/bash

#SBATCH --time=1-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=24
#SBATCH --mem=96Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=24

# Check if the correct number of arguments is provided
set -euo pipefail

#module load miniconda
#source activate chromsyn

FASTA="${1:?usage: $0 <FASTA>}"
TARGET=$(basename ${FASTA} .chr.fa)
FILE=$(realpath ${FASTA})
WD=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/busco_synteny
GENOMES=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/softmasked

echo "Working on ${TARGET}, file ${FASTA}"
export PYTHONWARNINGS="ignore::SyntaxWarning"

# Prep busco db 
BUSCO_DB=/project/coffea_pangenome/Software/Merondun/busco_downloads
LINEAGE=embryophyta_odb12
if [ -d "${BUSCO_DB}/lineages/${LINEAGE}" ]; then
        echo "BUSCO db ${LINEAGE} already present at ${BUSCO_DB}/lineages/${LINEAGE} – skipping download"
else
        busco --download ${LINEAGE} --download_path ${BUSCO_DB}
fi

mkdir -p ${WD}/work ${WD}/plotting_inputs
cd ${WD}/work

# Generate inputs
TELO_DIR=/project/coffea_pangenome/Software/Merondun/telociraptor/code
if [ -f ${TARGET}.chr.telomeres.tdt]; then
        echo "Telociraptor output exists for ${TARGET} – skipping"
else
        python ${TELO_DIR}/telociraptor.py seqin=${FILE} basefile=${FILE} i=-1 tweak=F telonull=T
fi

# busco 
if [ -f ${TARGET}.busco5.tsv]; then
        echo "BUSCO already ran on ${TARGET} - skipping"
else
        busco -f -o run_${TARGET} -i ${FILE} -l ${BUSCO_DB}/lineages/${LINEAGE} --cpu ${t} -m genome
        cp -v run_${TARGET}/run_${LINEAGE}/full_table.tsv ${TARGET}.busco5.tsv
        rm -rf run_${TARGET}*
fi

# repeat scores
if [ -f ${TARGET}.tidk.tsv]; then
        echo "TIDK already ran on ${TARGET} - skipping"
else
        tidk search --dir search --output ${TARGET} -s AACCCT ${FILE}
        cp -v search/${TARGET}_telomeric_repeat_windows.tsv ${TARGET}.tidk.tsv
fi

# Copy outputs
cp ${TARGET}.tidk.tsv ${TARGET}.chr.gaps.tdt ${TARGET}.busco5.tsv ${TARGET}.chr.telomeres.tdt ${TARGET}.chr.contigs.tdt ${WD}/plotting_inputs/
```

Merge the outputs and plot:

```bash
#!/bin/bash

#SBATCH --time=0-06:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

WD=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/busco_synteny
cd ${WD}

> busco.fofn > gaps.fofn > sequences.fofn > tidk.fofn

for i in $(cat Samples.list); do 
    echo -e "${i} ${WD}/plotting_inputs/${i}.busco5.tsv" >> busco.fofn
    echo -e "${i} ${WD}/plotting_inputs/${i}.chr.gaps.tdt" >> gaps.fofn
    echo -e "${i} ${WD}/plotting_inputs/${i}.chr.telomeres.tdt" >> sequences.fofn
    echo -e "${i} ${WD}/plotting_inputs/${i}.tidk.tsv" >> tidk.fofn
done 

Rscript ~/symlinks/software/chromsyn/chromsyn.R labelsize=1.5 opacity=0.4 pdfheight=8 pdfwidth=8
gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/printer -dNOPAUSE -dQUIET -dBATCH -sOutputFile=output.pdf chromsyn.pdf
```



## Invert Chr07 on HART001

Initial dotplots across all 14 accessions against sample HART001 identifies a consistent inversion on Chr07, which after closer inspection of the HiC heat map for HART001, is a putative misassembly. For consistency in this paper, I will reverse complement that contig (confirmed it is a long arm gapped with >50 Ns), so that we have a consistent Chr07 direction. 

![HART001_misassembly](/imgs/HART001_misassembly_inversion.png)

* isolate Chr07 second contig around big N-gap (≥51 Ns) → split into part1 / gap / part2 

* reverse‑complement part2 → rebuild Chr07 with inverted arm (p1 + gap + p2_rc) 

- rename rebuilt contig and verify lengths w/ seqkit 

- replace original Chr07 in HART001.fa with new inverted version → write HART001.Chr07sub.fa 

* Rearrange second contig of Chr07 in HART001. 

```bash
#split 
perl -0777 -pe '
s/^>.*\n//m;        
s/\s+//g;          
if (!/(N{51,})/) { die "ERROR: no N-run >=51 found\n"; }
my $gap = $1;
my ($a,$b) = split(/N{51,}/, $_, 2);
open(F1, ">part1.fa") or die $!;
print F1 ">Chr07_part1\n$a\n";
open(FG, ">gap.fa") or die $!;
print FG ">Chr07_gap\n$gap\n";
open(F2, ">part2.fa") or die $!;
print F2 ">Chr07_part2\n$b\n";
$_="";
' Chr07.fa


seqkit seq -r -p part2.fa > part2.rc.fa
p1=$(perl -ne 'next if /^>/; s/\s+//g; print' part1.fa)
gap=$(perl -ne 'next if /^>/; s/\s+//g; print' gap.fa)
p2rc=$(perl -ne 'next if /^>/; s/\s+//g; print' part2.rc.fa)

{
  echo ">Chr07_invB"
  printf "%s%s%s\n" "$p1" "$gap" "$p2rc"
} > Chr07_invB.fa

#verify
seqkit fx2tab -n -l Chr07.fa Chr07_invB.fa


perl -pe 'if($.==1){s/^>.*/>Chr07/}' Chr07_invB.fa > Chr07.repl.fa

perl -0777 -ne '
BEGIN {
  $repl = do {
    open my $fh, "<", "Chr07.repl.fa" or die $!;
    local $/; <$fh>
  };
}

# split into fasta records
@recs = split(/(?=>)/, $_);
for $r (@recs) {
  next if $r =~ /^\s*$/;

  if ($r =~ /^>(\S+)/) {
    $id = $1;
    if ($id eq "Chr07") {
      print $repl;
      $repl = "";  
    } else {
      print $r;
    }
  }
}
' HART001.fa > HART001.Chr07sub.fa

seqkit fx2tab -n -l HART001.fa | wc -l
seqkit fx2tab -n -l HART001.Chr07sub.fa | wc -l

seqkit fx2tab -n -l HART001.Chr07sub.fa | grep -w '^Chr07'
seqkit fx2tab -n -l Chr07.repl.fa
```

## HART069 "Butterfly"

Sample HART069 looks like it has a truncated Chr09 and an extra large Chr022, which appears as a misassembly in the dotplot against HART001. However, the juicer contact map for this sample does not show a straight story: instead it shows a potential chromosome-arm duplication, as it seems that the portion of chr09 missing ALSO matches chr22. 

This could require a phased genome approach to resolve, but archiving this to remember any issues with chr09 or chr22 in the future. 

![HART069 map](/imgs/HART069_inversion_chr22_chr09.png)

