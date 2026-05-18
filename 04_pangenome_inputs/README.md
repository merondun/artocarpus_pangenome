# Prepare Pangenome Inputs

For the pangenome, we will first start at the chromosome-level and build a collapsed pangenome for **chrN** using the 14 collapsed assemblies. 

First, rename the chrs so that they include sample and chromosome IDs, following the [panSN spec](https://github.com/pangenome/PanSN-spec). We will add '1' for each haplotype ID since these are collapsed primary assemblies. 

for `Chr01`, looks like:

```
zcat Chr01.fa.gz | grep '>' | head 
>HART001#1#Chr01
>HART030#1#Chr01
>HART050#1#Chr01
>HART069#1#Chr01
>HART032#1#Chr01
>HART033#1#Chr01
```

Submit:

```bash
#!/bin/bash
#SBATCH --time=0-05:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

set -euo pipefail

indir=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked
outdir=/project/coffea_pangenome/Breadfruit_Pangenome/01_assemblies

mkdir -p $outdir

for chrnum in $(seq 1 28); do
    chr=$(printf "Chr%02d" "$chrnum")
    outfile="$outdir/${chr}.fa"
    : > "$outfile"

    for sample in $(cat ../Samples.list); do 

        awk -v chr="$chr" -v sample="$sample" '
            /^>/ {
                header=$0
                sub(/^>/, "", header)

                split(header, fields, /[[:space:]]+/)

                if (fields[1] == chr) {
                    print ">" sample "#" 1 "#" chr
                    keep=1
                } else {
                    keep=0
                }
                next
            }
            keep { print }
        ' ${indir}/${sample}.fa >> "$outfile"
    done

    bgzip -f $outfile;
    samtools faidx "$outfile.gz"
    nchrs=$(zcat ${outfile} | grep '>' | wc -l)
    echo "there are ${nchrs} chromosomes for ${chr}"
done

```

There should be 14 representatives for each chr:

```
cat slurm-20878430.out 
there are 14 chromosomes for Chr01
there are 14 chromosomes for Chr02
there are 14 chromosomes for Chr03
there are 14 chromosomes for Chr04
there are 14 chromosomes for Chr05
there are 14 chromosomes for Chr06
there are 14 chromosomes for Chr07
there are 14 chromosomes for Chr08
there are 14 chromosomes for Chr09
there are 14 chromosomes for Chr10
there are 14 chromosomes for Chr11
there are 14 chromosomes for Chr12
there are 14 chromosomes for Chr13
there are 14 chromosomes for Chr14
there are 14 chromosomes for Chr15
there are 14 chromosomes for Chr16
there are 14 chromosomes for Chr17
there are 14 chromosomes for Chr18
there are 14 chromosomes for Chr19
there are 14 chromosomes for Chr20
there are 14 chromosomes for Chr21
there are 14 chromosomes for Chr22
there are 14 chromosomes for Chr23
there are 14 chromosomes for Chr24
there are 14 chromosomes for Chr25
there are 14 chromosomes for Chr26
there are 14 chromosomes for Chr27
there are 14 chromosomes for Chr28
```
