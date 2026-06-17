# 07_annotation/

Annotate using egapx with all available breadfruit data: 6 long read ISOseq libraries from HART032 (n=4) and 2 from H6 (n=2). As well as 16 short read libraries from the SRA.

```
cat Isoreads.txt 
H6ML    /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/isoseq/H6__ML.fasta
H6YL    /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/isoseq/H6__YL.fasta
HFR     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/isoseq/HART032__FR.fasta
HMF     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/isoseq/HART032__MF.fasta
HML     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/isoseq/HART032__ML.fasta
HYL     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/isoseq/HART032__YL.fasta

cat Shortreads.txt 
S1      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997516_1.fastq
S1      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997516_2.fastq
S2      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997517_1.fastq
S2      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997517_2.fastq
S3      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997518_1.fastq
S3      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997518_2.fastq
S4      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997519_1.fastq
S4      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997519_2.fastq
S5      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997520_1.fastq
S5      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997520_2.fastq
S6      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997521_1.fastq
S6      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997521_2.fastq
S7      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997522_1.fastq
S7      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997522_2.fastq
S8      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997523_1.fastq
S8      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997523_2.fastq
S9      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997524_1.fastq
S9      /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997524_2.fastq
S10     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997525_1.fastq
S10     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997525_2.fastq
S11     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997526_1.fastq
S11     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997526_2.fastq
S12     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997527_1.fastq
S12     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997527_2.fastq
S13     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997536_1.fastq
S13     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997536_2.fastq
S14     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997537_1.fastq
S14     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997537_2.fastq
S15     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997538_1.fastq
S15     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997538_2.fastq
S16     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997539_1.fastq
S16     /90daydata/coffea_pangenome/scratch/annotation_egapx_pangenome/rnaseq/external_SRR5997539_2.fastq
```

These are the accessions for the publicly available shortread data that I will use for annotation:

| Run        | Bases      | geo_loc_name | Library Name | Organism                                     | Cultivar      | tissue   | dev_stage |
| ---------- | ---------- | ------------ | ------------ | -------------------------------------------- | ------------- | -------- | --------- |
| SRR5997516 | 3319449840 | USA: Hawaii  | RNA26        | Artocarpus altilis                           | Meisaip       | Leaf     | Mature    |
| SRR5997517 | 3465803486 | USA: Hawaii  | RNA25        | Artocarpus altilis                           | Meisei        | Leaf     | Mature    |
| SRR5997518 | 3287302954 | USA: Hawaii  | RNA17        | Artocarpus altilis x Artocarpus mariannensis | Unk (Huehue)  | Leaf     | Mature    |
| SRR5997519 | 1369060858 | USA: Hawaii  | RNA40        | Artocarpus altilis x Artocarpus mariannensis | Faine         | Leaf     | Mature    |
| SRR5997520 | 1797107342 | USA: Hawaii  | RNA36        | Artocarpus altilis x Artocarpus mariannensis | Ulu afa elise | Perianth | Mature    |
| SRR5997521 | 3894145294 | USA: Hawaii  | RNA32        | Artocarpus altilis x Artocarpus mariannensis | Ulu afa       | Perianth | Mature    |
| SRR5997522 | 4491655032 | USA: Hawaii  | EW1          | Artocarpus altilis                           | Ulu fiti      | Perianth | Mature    |
| SRR5997523 | 3617708900 | USA: Hawaii  | RNA10        | Artocarpus altilis x Artocarpus mariannensis | Rotuma        | Leaf     | Mature    |
| SRR5997524 | 4021107344 | USA: Hawaii  | RNA24        | Artocarpus altilis                           | Samoan 2      | Leaf     | Mature    |
| SRR5997525 | 4842006862 | USA: Hawaii  | RNA21        | Artocarpus altilis                           | Karawa        | Perianth | Mature    |
| SRR5997526 | 4739479742 | USA: Hawaii  | RNA16        | Artocarpus altilis x Artocarpus mariannensis | Midolab       | Leaf     | Mature    |
| SRR5997527 | 2246774694 | USA: Hawaii  | RNA48        | Artocarpus altilis x Artocarpus mariannensis | Luthar        | Leaf     | Mature    |
| SRR5997536 | 4718087134 | USA: Hawaii  | RNA38        | Artocarpus altilis                           | Kea           | Leaf     | Mature    |
| SRR5997537 | 4970992548 | USA: Hawaii  | EW3          | Artocarpus altilis                           | Toneno        | Perianth | Mature    |
| SRR5997538 | 4203056016 | USA: Hawaii  | RNA39        | Artocarpus altilis                           | Puupuu        | Leaf     | Mature    |
| SRR5997539 | 3132807698 | USA: Hawaii  | RNA49        | Artocarpus altilis                           | Rotuma        | Leaf     | Mature    |

Then submit this with e.g.: `sbatch 01_Annotation_eGAPx.sh HART050` 

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

module load egapx/0.5.2

# Variables for genome to annotate, and which isoseq reads to use 
ID=${1:?Provide ID as first argument}

# Static
RUN="${ID}_full"
WD=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/annotation
GENOMES=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked
echo "WORKING ON ${RUN}"

# Create yaml for run 
YAML=${WD}/${RUN}.yaml
cat > "$YAML" <<EOF
genome: ${GENOMES}/${ID}.fa
taxid: 194251
long_reads: ${WD}/Isoreads.txt
short_reads: ${WD}/Shortreads.txt
cmsearch:
  enabled: false
trnascan:
  enabled: false
EOF

# Run 
egapx.py ${YAML} -e slurm -w ${WD}/work/${RUN} -o ${WD}/out/${RUN}
```

Afterwards, go to the output folder and extract the gtfs, just count genes for sanity first:

```bash
mkdir -p ../gtfs
for SAMPLE in */; do ID=$(echo $SAMPLE | sed 's@_full/@@g'); cp ${SAMPLE}/complete.genomic.gtf ../gtfs/${ID}.gtf; ngenes=$(awk '$3 == "gene"' ../gtfs/${ID}.gtf | wc -l); echo -e "${ID}\t${ngenes}"; done
```

Genes:

```bash
H6      30991
HART001 31763
HART030 31523
HART032 31440
HART033 31357
HART038 31443
HART046 31098
HART049 30497
HART050 31583
HART053 31128
HART069 31334
ZZ3     30824
ZZ7     30750
ZZ9     31008
```

