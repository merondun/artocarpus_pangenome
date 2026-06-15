# 03_repeats/

* Repeats are annotated with EarlGrey on all assemblies. This then gathers summary files and merges outputs into combined repeat family, summary, and divergence tables.
* Processes merged to aggregate repeat classes and generate high‑level repeat composition plots and summary tables for all accessions.

Output:



![repeats](/imgs/20260612_RepeatsHighLevelSummary.png)



___

Run **earlgrey v7.2.2**:

```bash
#!/bin/bash

#SBATCH --time=14-00:00:00   
#SBATCH --cpus-per-task=16
#SBATCH --mem=384Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=16

#module load miniconda
#source activate earlgrey722

SAMPLE=$1
WD=/90daydata/coffea_pangenome/scratch/repeats
ASM=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/primary_asm
OUTDIR=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/repeats

cd ${WD}

FASTA="${ASM}/${SAMPLE}.fa" 

echo -e "\e[43m~~~~ Starting repeat annotation for ${SAMPLE} ~~~~\e[0m"
# Run earlgrey with eudicotyledons repeatmasker search time, generating soft-masked genome and run helitrons. 
earlGrey -r eudicotyledons -d yes -t ${t} -g ${FASTA} -q yes -s ${SAMPLE} -o ${OUTDIR}/
```

And copy:

```bash
DIR=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/repeats
REP=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/repeats/output

for SAMPLE in $(cat Pangenome.list) ; do 
cp ${DIR}/${SAMPLE}_EarlGrey/${SAMPLE}_summaryFiles/* ${REP}/
done 
```

Add accession to each output:

```bash
for SAMPLE in $(ls *.familyLevelCount.txt | sed 's/.familyLevelCount.txt//g'); do 
echo "${SAMPLE}"
awk -v s=${SAMPLE} '{OFS="\t"}{print $0, s}' ${SAMPLE}.familyLevelCount.txt > ${SAMPLE}.families.out
awk -v s=${SAMPLE} '{OFS="\t"}{print $0, s}' ${SAMPLE}.highLevelCount.txt > ${SAMPLE}.summary.out
awk -v s=${SAMPLE} '{OFS="\t"}{print $0, s}' ${SAMPLE}_divergence_summary_table.tsv > ${SAMPLE}.divergence.out
done 

mergem *families.out > Repeat_Families.txt
mergem *summary.out > Repeat_Summaries.txt
mergem *divergence.out > Divergence_Summaries.txt
```

### Plot

```R
setwd('/project/coffea_pangenome/Artocarpus/Pangenome_Paper/repeats/output')
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(meRo) #devtools::install_github('merondun/meRo')
library(vegan)
library(broom)
library(ggrepel)
library(caper)

# Add metadata information
md <- read_tsv('~/artocarpus_pangenome/samples.info')

##### High Level #####
high_level <- read_tsv('~/artocarpus_pangenome/03_repeats/Repeat_Summaries.txt') 
names(high_level) <- c('Classification','Coverage','Count','Proportion','Gen','Distinct_Classifications','Sample')
fl <- left_join(high_level,md)
fl$Sample <- factor(fl$Sample,levels=md$Sample)
fl <- fl %>%
  arrange(Sample) %>%
  mutate(Group = factor(Group, levels = unique(Group)))

# aggregate nested
fl <- fl %>%
  mutate(ClassSimple = gsub("-nested", "", Classification))
fl_sum <- fl %>%
  group_by(Sample, Group, ClassSimple) %>%
  summarise(Proportion = sum(Proportion), .groups="drop") %>% 
  filter(ClassSimple != "Total Interspersed Repeat")

fl_sum$ClassSimple <- factor(fl_sum$ClassSimple,levels=c('Non-Repeat','Unclassified','Other (Simple Repeat, Microsatellite, RNA)','DNA','Penelope','Rolling Circle','LTR','LINE','SINE'))
cols <- fl_sum %>% dplyr::select(ClassSimple) %>% distinct %>% mutate(Color = brewer.pal(9,'Paired'))

# ltr labs
fl_labels <- fl_sum %>%
  filter(ClassSimple == "LTR") %>%
  mutate(
    label = paste0(round(Proportion, 1), "%"),
    text_color = ifelse(Proportion > 8, "black", "black")
  )


# Plot landscape 
all <- fl_sum %>% 
  #mutate(Coverage = Coverage / 1e6) %>% 
  pivot_longer(c(Proportion)) %>%
  #filter( !(name == 'Distinct_Classifications' & (Classification == 'Unclassified' | Classification == "Other (Simple Repeat, Microsatellite, RNA)")) ) %>% 
  ggplot(aes(y = Sample, x = value, fill = ClassSimple)) +
  geom_bar(stat = 'identity', position = position_stack()) +
  #add LTR percent labels
  geom_text(
    data = fl_labels,
    aes(y = Sample, x = Proportion, label = label),
    position = position_stack(vjust = 0.5),
    color = fl_labels$text_color,
    size = 2.5
  ) +
  theme_bw() +
  facet_grid(Group ~ name, scales = 'free', space = 'free_y') +
  scale_fill_manual(values = cols$Color, breaks = cols$ClassSimple) +
  theme(strip.text.y = element_text(angle = 0)) +
  ylab('') + xlab('Distinct Classifications') +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

all
ggsave('~/artocarpus_pangenome/03_repeats/20260612_RepeatsHighLevelSummary.pdf',
       all,dpi=300,height=5,width=6.5)

fl_sum %>% dplyr::select(Sample,Group,ClassSimple,Proportion) %>% pivot_wider(names_from = ClassSimple,values_from=Proportion)
# Sample  Group     DNA  LINE   LTR `Non-Repeat` Other (Simple Repeat…¹ Penelope `Rolling Circle`    SINE Unclassified
# <fct>   <fct>   <dbl> <dbl> <dbl>        <dbl>                  <dbl>    <dbl>            <dbl>   <dbl>        <dbl>
#   1 HART001 A. alt…  4.05 0.387  42.0         32.0                   2.33        0            0.415 0.00151         19.8
# 2 HART030 A. alt…  3.87 0.402  40.2         32.1                   2.31        0            0.896 0.00116         20.9
# 3 HART050 A. alt…  3.88 0.325  42.4         31.4                   2.18        0            0.426 0.00174         20.2
# 4 HART069 A. alt…  4.15 0.346  42.4         32.0                   2.34        0            0.411 0.00349         19.1
# 5 HART032 A. alt…  3.92 0.412  40.9         32.0                   2.03        0            0.783 0.00613         20.7
# 6 HART033 A. alt…  3.64 0.303  40.0         33.0                   2.66        0            0.464 0.00336         21.1
# 7 HART038 A. alt…  4.11 0.253  40.7         31.9                   2.24        0            0.411 0.00216         21.4
# 8 ZZ3     A. alt…  4.96 0.502  37.6         33.5                   2.19        0            0.381 0.00156         21.5
# 9 ZZ7     A. alt…  4.57 0.390  39.2         33.6                   2.09        0            0.407 0.00146         20.6
# 10 ZZ9     A. alt…  4.20 0.419  38.5         33.6                   2.31        0            0.403 0.00605         21.2
# 11 HART046 A. alt…  4.25 0.441  38.7         33.1                   3.11        0            0.685 0.00174         20.6
# 12 HART049 A. alt…  4.39 0.411  39.3         32.8                   2.69        0            0.539 0.00543         20.7
# 13 HART053 A. alt…  4.08 0.309  39.5         33.0                   2.42        0            0.570 0.00158         20.8
# 14 H6      A. alt…  4.56 0.376  40.2         33.0                   2.33        0            0.412 0.00220         20.0

write.csv(fl_sum %>% dplyr::select(Sample,Group,ClassSimple,Proportion),'~/artocarpus_pangenome/03_repeats/Repeat_proportions_coverage_summarized.csv',quote = F,row.names = F)

```



