# Prepare Pangenome Inputs

For the pangenome, we will first start at the chromosome-level and build a collapsed pangenome for **chrN** using the 14 collapsed assemblies. 

First, rename the chrs so that they include sample and chromosome IDs, following the [panSN spec](https://github.com/pangenome/PanSN-spec) for non-haplotype phased:

for `Chr01`, looks like:

```
zcat Chr01.fa.gz | grep '>' | head
>H6#Chr01
>HART001#Chr01
>HART030#Chr01
>HART032#Chr01
>HART033#Chr01
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

mkdir -p "$outdir"

for chr in $(seq -w 1 28); do
    outfile="$outdir/Chr${chr}.fa"
    : > "$outfile"

    for fa in "$indir"/*.chr.fa; do
        sample=$(basename "$fa" .chr.fa)

        awk -v chr="Chr${chr}" -v sample="$sample" '
            /^>/ {
                header=$0
                sub(/^>/, "", header)
                if (header ~ chr) {
                    print ">" sample "#" chr
                    keep=1
                } else {
                    keep=0
                }
                next
            }
            keep { print }
        ' "$fa" >> "$outfile"
    done

    gzip -f "$outfile"
done
```

