# 01_qa_qc_genomescope/

* basic QC: seqkit stats per sample → merge into Read_Lengths.txt 

* run genomescope2 + smudgeplot pipeline (FastK → Histex → GS2 → smudgeplot) 
* extract GS2 metrics (het, haploid size, unique content) into .info files 



Outputs:

Assembly inputs:

| Sample  | Ploidy | Group                | HiFi   | HiFiLength | HiC      | Name          | Color   | HapCoverage |
| ------- | ------ | -------------------- | ------ | ---------- | -------- | ------------- | ------- | ----------- |
| HART001 | 2      | A. altilis 2N        | 59.32  | 11206.3    | 50.46    | Maafala       | #0072B2 | 29.66       |
| HART030 | 2      | A. altilis 2N        | 96.64  | 18338.2    | 17.63    | Huero         | #0072B2 | 48.32       |
| HART050 | 2      | A. altilis 2N        | 38.96  | 17753.1    | 18.07    | Kukumu tasi   | #0072B2 | 19.48       |
| HART069 | 2      | A. altilis 2N        | 72.95  | 12674.4    | 64.45    | Ulu Fiti      | #0072B2 | 36.475      |
| HART032 | 3      | A. altilis 3N        | 88.25  | 15382.2    | 26.43    | Hamoa         | #56B4E9 | 29.41667    |
| HART033 | 3      | A. altilis 3N        | 139.76 | 14764.85   | 36.96617 | Patara        | #56B4E9 | 46.58667    |
| HART038 | 3      | A. altilis 3N        | 42.63  | 18956.6    | 29.35    | Lemai         | #56B4E9 | 14.21       |
| ZZ3     | 2      | A. altilis hybrid 2N | 29.22  | 18949.4    | 0        | Ulu hamoa     | #D55E00 | 14.61       |
| ZZ7     | 2      | A. altilis hybrid 2N | 28.6   | 20738.2    | 0        | Ulu afa elise | #D55E00 | 14.3        |
| ZZ9     | 2      | A. altilis hybrid 2N | 30.02  | 19115.8    | 0        | Ulu afa       | #D55E00 | 15.01       |
| HART046 | 3      | A. altilis hybrid 3N | 111.19 | 14132.65   | 37.7     | Midolab       | #E69F00 | 37.06333    |
| HART049 | 3      | A. altilis hybrid 3N | 92.85  | 8248.15    | 3.47     | Meinpohnsakar | #E69F00 | 30.95       |
| HART053 | 3      | A. altilis hybrid 3N | 78.23  | 15006      | 32.63    | Nahnmwal      | #E69F00 | 26.07667    |
| H6      | 4      | A. altilis hybrid 4N | 87.41  | 15742.8    | 35.28212 | Meikole       | #009E73 | 21.8525     |

And both genomescope2 kmer spectra and ploidy smudgeplots:

![smudgeplots](/imgs/smudge_genomescope.png)



___



Read lengths:

```bash
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --partition=short

WD=/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies

SAMPLE=$1

seqkit stats -a ${WD}/${SAMPLE}/${SAMPLE}.HiFi.fastq.gz > ${SAMPLE}.stats
```

Afterwards:

```bash
mergem *stats | sed 's@.*\/@@g; s/.HiFi.fastq.gz//g; s/[,()]//g' > Read_Lengths.txt
```

And heterozygosity from genomescope2 (below):

```bash
for i in $(ls *summary.txt | sed 's/_summary.txt//g'); do grep 'Het' ${i}_summary.txt | sed 's/.*)//g' | sed 's/%.*//g' | awk -v i=${i} '{OFS="\t"}{print $0, i}' > ${i}.het; done

cat *het | sed 's/ //g' > Heterozygosity.txt
```

And Stats if needed:

```bash
#!/bin/bash
#SBATCH --time=24:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem=30Gb
#SBATCH --partition=short

FILE="$1"

# Get stats on bases / reads, will output to slurm 
bbduk.sh in=${FILE}

# Extract 
BASES=$(grep 'Input:' slurm-$SLURM_JOB_ID.out | awk '{print $4}')
READS=$(grep 'Input:' slurm-$SLURM_JOB_ID.out | awk '{print $2}')
MD5SUM=$(md5sum ${FILE})
echo -e "$MD5SUM\t$BASES\t$READS" > ${FILE}.info
```

## K-mer Based Genome Size Estimation: Genomescope

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --cpus-per-task=24
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome
 
set -euo pipefail

module load miniconda
source activate puzzler200 

SAMPLE=${1:?Provide SAMPLE id as first argument}

# Read tsv
IFS=$'\t,' read -r _ PLOIDY HIFI < <(
    awk -F'[\t,]' -v sample="$SAMPLE" '$1 == sample {print $0}' smudge_info.tsv
)

# Where our analyses are run 
WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/smudges
echo -e "=======================================================================\nParameters for sample: ${SAMPLE} \nPLOIDY: ${PLOIDY}\nHIFI: ${HIFI}\n=======================================================================\n"

# Create directory to save all results 
mkdir -p ${WD}/results_${SAMPLE}
cd ${WD}/results_${SAMPLE}
echo "Working on ${SAMPLE} for file ${HIFI}" 

# This is a kmer value, in general 31 works well for a range of eukaryotes
K=31

# First, generate your kmer table, this is memory and time intensive 
FastK -v -t1 -k${K} -M70 -T24 -NFastK_Table_k${K}_${SAMPLE} ${HIFI}

# Now, generate a coverate histogram from that table 
Histex -G FastK_Table_k${K}_${SAMPLE} > ${SAMPLE}_k${K}.hist

# Run genomescope2 on that histogram table 
genomescope2 --input ${SAMPLE}_k${K}.hist --output . --ploidy ${PLOIDY} --kmer_length ${K} --name_prefix ${SAMPLE}_k${K}

conda deactivate
source activate smudge 

# Create directory to save all results
cd ${WD}/results_${SAMPLE}

# Find all k-mer pairs in the dataset using hetmer module
smudgeplot hetmers -L 12 -t 4 -tmp . -o ${SAMPLE}_kmerpairs --verbose FastK_Table_k${K}_${SAMPLE}

# this now generated `data/Scer/kmerpairs_text.smu` file;
# it's a flat file with three columns; covB, covA and freq (the number of k-mer pairs with these respective coverages)

# use the .smu file to infer ploidy and create smudgeplot
smudgeplot all --format pdf -o ${SAMPLE}_smudge ${SAMPLE}_kmerpairs.smu

# Afterwards, you will want to clean up that massive file 
#rm FastK*

```

Extract info from genomescope: 

```bash
### Extract Stats
for i in $(ls *summary.txt | sed 's/_k31_summary.txt//g' ); do 

echo "Processing ${i}"

HETMIN=$(grep 'Het' ${i}*_summary.txt | sed 's/.*)//g' | awk '{print $1}') 
HETMAX=$(grep 'Het' ${i}*_summary.txt | sed 's/.*)//g' | awk '{print $2}') 
HAPMIN=$(grep 'Hap' ${i}*_summary.txt | sed 's/.*)//g' | awk '{print $4}' | sed 's/,//g')
HAPMAX=$(grep 'Hap' ${i}*_summary.txt | sed 's/.*)//g' | awk '{print $6}' | sed 's/,//g')
UNIMIN=$(grep 'Unique' ${i}*_summary.txt | sed 's/.*)//g' | awk '{print $4}' | sed 's/,//g')
UNIMAX=$(grep 'Unique' ${i}*_summary.txt | sed 's/.*)//g' | awk '{print $6}' | sed 's/,//g')

echo -e "Sample\tGS_HetMin\tGS_HetMax\tGS_SizeMin\tGS_SizeMax\tGS_UniqueMin\tGS_UniqueMax" > ${i}.info
echo -e "${i}\t${HETMIN}\t${HETMAX}\t${HAPMIN}\t${HAPMAX}\t${UNIMIN}\t${UNIMAX}" >> ${i}.info

done 
```

Plot:

```R
setwd('/project/coffea_pangenome/Artocarpus/Pangenome_Paper/smudges/outputs/')
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
library(viridis)
library(ggtext)

md <- read.table('../../samples.info',sep='\t',header = TRUE,comment.char = '') %>% as_tibble
grpcol <- md %>% distinct(Group, Color) %>% deframe()

info_dat <- NULL
hist_dat <- NULL
files <- list.files('.',pattern = '*_k31.info')
for (file in files) {
  id = gsub('_k31.info','',file)
  cat('Processing: ',id,'\n')
  info <- read_tsv(file)
  info_dat <- rbind(info_dat,info)
  hist <- read_tsv(paste0(id,'_k31.hist'),col_names=F)
  hist <- hist %>% mutate(Sample = id)
  hist_dat <- rbind(hist_dat,hist)
}

names(hist_dat) <- c('Bin','Coverage','Sample')
hm <- left_join(hist_dat,md)
hm$Sample <- factor(hm$Sample,levels=unique(md$Sample))

# limits
xlims <- hm %>% group_by(Sample,Group) %>% summarize(xmin=5,xmax=HapCoverage*Ploidy+HapCoverage*2.5) %>% distinct

hm2 <- hm %>%
  left_join(xlims, by = c("Sample", "Group")) %>%
  mutate(ID = paste0(Sample, " (<i>", Group, "</i>)"),
         ID = factor(ID, levels = paste0(levels(Sample), " (<i>", Group[match(levels(Sample), Sample)], "</i>)"))) %>%
  filter(Bin >= xmin, Bin <= xmax)

hist_plot <- hm2 %>%
  ggplot(aes(x = Bin, y = Coverage, fill = Group)) +
  geom_bar(stat='identity') +
  facet_wrap(~ ID, scales = "free", nrow = 7, ncol = 2) +
  scale_fill_manual(values = grpcol) +
  theme_bw(base_size=9) +
  labs(x='Coverage Peak',y='')+
  theme(
    strip.text = element_markdown(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
hist_plot

ggsave('~/symlinks/pan/figures/20260504_Genomescope_Histograms.pdf',hist_plot,height=8,width=4.5,dpi=300)

```

