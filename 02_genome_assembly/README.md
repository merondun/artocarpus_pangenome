# 02_genome_assembly/

Primary collapsed assembly proceeds with [puzzler v2.0.0](https://github.com/merondun/puzzler) from HiFi and HiC. 

Outputs: assemblies & statistics:

| Accession | Group                | HiFi Effective  Coverage | Size  | CN50 | SN50 | Completeness | Chr Anchored | BUSCO_Complete | BUSCO_singlecopy |
| --------- | -------------------- | ------------------------ | ----- | ---- | ---- | ------------ | ------------ | -------------- | ---------------- |
| HART001   | A. altilis 2N        | 29.7                     | 843.3 | 15.4 | 30.9 | 87.5         | 99.8         | 98.3           | 69.2             |
| HART030   | A. altilis 2N        | 48.3                     | 836.5 | 12.7 | 29.8 | 86.8         | 99.4         | 98.6           | 60.4             |
| HART050   | A. altilis 2N        | 19.5                     | 855.5 | 15.8 | 31.4 | 86.6         | 99.9         | 98.2           | 59.1             |
| HART069   | A. altilis 2N        | 36.5                     | 838.3 | 13.8 | 30.3 | 86.8         | 99.7         | 98.6           | 61               |
| HART032   | A. altilis 3N        | 29.4                     | 837.5 | 19.8 | 30.3 | 84.5         | 99.9         | 98.2           | 59.7             |
| HART033   | A. altilis 3N        | 46.6                     | 833.7 | 19.7 | 29.9 | 82           | 99.7         | 98.6           | 60.7             |
| HART038   | A. altilis 3N        | 14.2                     | 838.2 | 15.4 | 30.2 | 84.2         | 99.7         | 98.6           | 60.6             |
| ZZ3       | A. altilis hybrid 2N | 14.6                     | 783.6 | 18.3 | 28.6 | 91.5         | 99.8         | 98.7           | 61.7             |
| ZZ7       | A. altilis hybrid 2N | 14.3                     | 778.2 | 18   | 28.4 | 94.1         | 99.9         | 98.7           | 61.3             |
| ZZ9       | A. altilis hybrid 2N | 15                       | 784.5 | 21.2 | 28.8 | 92.3         | 99.6         | 98.6           | 61.9             |
| HART046   | A. altilis hybrid 3N | 37.1                     | 800.4 | 19.3 | 28.6 | 70.7         | 99.5         | 98.5           | 61.7             |
| HART049   | A. altilis hybrid 3N | 31                       | 790   | 14.4 | 28.7 | 72.9         | 99.6         | 98.5           | 62.7             |
| HART053   | A. altilis hybrid 3N | 26.1                     | 796.6 | 18.3 | 29.1 | 73.6         | 99.7         | 98.6           | 61.6             |
| H6        | A. altilis hybrid 4N | 21.9                     | 798.2 | 18.3 | 28.5 | 70.2         | 99.9         | 98.5           | 61.8             |



Final HiC contact maps, pre- and post- juicer contacts, and merqury completeness and kmer spectra:

![HiC and merqury](/imgs/20260501_contacts_compiled.png) 



Percent of assembly present in large scaffolds in an Nx plot:

![Nx plot](/imgs/20260615_Nx_Plot.png)



___



Most samples (particularly diploids), were ran end-to-end without any issues - including a manual curation step:

```bash
#!/bin/bash

#SBATCH --time=8-00:00:00   
#SBATCH --cpus-per-task=48
#SBATCH --mem=512Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=48

#module load miniconda
#source activate puzzler200

SAMPLE=${1:?Provide SAMPLE id as first argument}

puzzler -s ${SAMPLE} -m samples.tsv --threads 48 --mem 508
```

From this `samples.tsv`:

| sample  | runtime | container | wd                                                           | hifi                                                         | hic_r1                                                       | hic_r2                                                       | num_chrs | reference                                                    | hom_cov | fcs_db                                          | fcs_taxid | busco_lineage     | busco_database                                              |
| ------- | ------- | --------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | ------- | ----------------------------------------------- | --------- | ----------------- | ----------------------------------------------------------- |
| HART030 | conda   | NA        | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART030.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART030.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART030.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/HART001.chr.fa | 108     | /project/coffea_pangenome/Software/Merondun/fcs | 194251    | embryophyta_odb12 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART050 | conda   | NA        | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART050.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART050.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART050.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/HART001.chr.fa | 42      | /project/coffea_pangenome/Software/Merondun/fcs | 194251    | embryophyta_odb12 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART069 | conda   | NA        | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART069.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART069.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART069.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/HART001.chr.fa | 85      | /project/coffea_pangenome/Software/Merondun/fcs | 194251    | embryophyta_odb12 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART032 | conda   | NA        | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART032.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART032.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART032.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/HART001.chr.fa | 102     | /project/coffea_pangenome/Software/Merondun/fcs | 194251    | embryophyta_odb12 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART033 | conda   | NA        | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART033.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART033.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART033.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/HART001.chr.fa | 162     | /project/coffea_pangenome/Software/Merondun/fcs | 194251    | embryophyta_odb12 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART038 | conda   | NA        | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART038.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART038.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART038.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/HART001.chr.fa | 48      | /project/coffea_pangenome/Software/Merondun/fcs | 194251    | embryophyta_odb12 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| ZZ3     | conda   | NA        | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ZZ3.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ZZ3.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ZZ3.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/HART001.chr.fa | 36      | /project/coffea_pangenome/Software/Merondun/fcs | 194251    | embryophyta_odb12 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| ZZ7     | conda   | NA        | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ZZ7.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ZZ7.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ZZ7.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/HART001.chr.fa | 36      | /project/coffea_pangenome/Software/Merondun/fcs | 194251    | embryophyta_odb12 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| ZZ9     | conda   | NA        | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ZZ9.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ZZ9.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ZZ9.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/HART001.chr.fa | 36      | /project/coffea_pangenome/Software/Merondun/fcs | 194251    | embryophyta_odb12 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART046 | conda   | NA        | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART046.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART046.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART046.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/HART001.chr.fa | 132     | /project/coffea_pangenome/Software/Merondun/fcs | 194251    | embryophyta_odb12 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART049 | conda   | NA        | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART049.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART049.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART049.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/HART001.chr.fa | 111     | /project/coffea_pangenome/Software/Merondun/fcs | 194251    | embryophyta_odb12 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART053 | conda   | NA        | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART053.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART053.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART053.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/HART001.chr.fa | 90      | /project/coffea_pangenome/Software/Merondun/fcs | 194251    | embryophyta_odb12 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| H6      | conda   | NA        | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/H6.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/H6.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/H6.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/HART001.chr.fa | 104     | /project/coffea_pangenome/Software/Merondun/fcs | 194251    | embryophyta_odb12 | /project/coffea_pangenome/Software/Merondun/busco_downloads |



## Merqury Histograms

- collect Merqury k‑mer copy‑number histograms for all accessions and combine into a unified dataset 

Just grab the merqury histos:

```
for i in $(cat Pangenome.list); do cp ${i}/08_merqury/merq.${i}.spectra-cn.hist merqury_spectra/;done
```

And plot:

```R
setwd('/project/coffea_pangenome/Artocarpus/Pangenome_Paper/merqury_spectra/')
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
library(viridis)
library(ggtext)

md <- read.table('../samples.info',sep='\t',header = TRUE,comment.char = '') %>% as_tibble
grpcol <- md %>% distinct(Group, Color) %>% deframe()

hist_dat <- NULL
files <- list.files('.',pattern = '*hist')
for (file in files) {
  id = gsub('merq.','',gsub('.spectra-cn.hist','',file))
  cat('Processing: ',id,'\n')
  hist <- read_tsv(file) %>% mutate(Sample = id)
  hist_dat <- rbind(hist_dat,hist)
}

hm <- left_join(hist_dat,md)
hm$Sample <- factor(hm$Sample,levels=unique(md$Sample))

# limits
xlims <- hm %>% group_by(Sample,Group) %>% summarize(xmin=5,xmax=HapCoverage*Ploidy+HapCoverage*2.5) %>% distinct

hm2 <- hm %>%
  left_join(xlims, by = c("Sample", "Group")) %>%
  mutate(ID = paste0(Sample, " (<i>", Group, "</i>; ",round(Merqury_Complete,1),'%)')) %>% 
  filter(kmer_multiplicity >= xmin, kmer_multiplicity <= xmax ) %>% 
  mutate(Copies = factor(Copies,levels = c("read-only", "1"     ,    "2"        , "3" ,        "4"    ,     ">4")))
ord <- hm2 %>% arrange(Sample) %>%  distinct(Sample,ID)
hm2$ID <- factor(hm2$ID,levels=ord$ID)

hist_plot <- hm2 %>%
  ggplot(aes(x = kmer_multiplicity, y = Count, fill = Copies)) +
  geom_bar(stat='identity',alpha=0.6) +
  facet_wrap(~ ID, scales = "free", nrow = 7, ncol = 2) +
  scale_fill_manual(values = c('grey40',viridis(5))) +
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

ggsave('~/symlinks/pan/figures/20260505_Merqury_Histograms.pdf',hist_plot,height=8,width=4.5,dpi=300)

```

## Plot Nx Curve

• Computes Nx curves for all genome assemblies, following the ideas from [here](https://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity) 

```R
setwd('/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked')
library(Biostrings)
library(ggplot2)
library(tidyverse)
library(purrr)

md <- read_tsv('~/artocarpus_pangenome/samples.info')

files <- list.files('.', pattern="\\.fa$", full.names=TRUE)
files <- files[!grepl('chr',files)]

# Compute Nx curve for a vector of contig lengths
compute_Nx_curve <- function(lengths) {
  lengths <- sort(lengths, decreasing = TRUE)
  total <- sum(lengths)
  cumcov <- cumsum(lengths) / total * 100
  
  xvals <- 0:100
  nx <- sapply(xvals, function(x) {
    idx <- which(cumcov >= x)[1]
    if (is.na(idx)) return(0)
    lengths[idx]
  })
  
  data.frame(x = xvals, Nx = nx)
}

# Load contig lengths
get_contig_lengths <- function(fasta_file) {
  seqs <- readDNAStringSet(fasta_file)
  width(seqs)
}

# Loop through files and compute Nx for each
nx_list <- list()

for (i in seq_along(files)) {
  f <- files[i]
  contig_lengths <- get_contig_lengths(f)
  nx_df <- compute_Nx_curve(contig_lengths)
  
  # Add identifiers
  nx_df$Assembly <- basename(f)
  
  nx_list[[i]] <- nx_df
}

# Combine into a single data frame
all_nx <- bind_rows(nx_list) %>% 
  mutate(Assembly = gsub('.fa','',Assembly)) %>% 
  left_join(.,md %>% dplyr::select(Assembly=Sample,Group,Color))

# Plot Nx curves on one figure
nx <- ggplot(all_nx, aes(x = x, y = Nx, group = Assembly, color = Group)) +
  geom_line(size = 0.25) +
  theme_bw(base_size=11) +
  xlab("Percent Assembly Covered") +
  ylab("Nx (Contig Length)") +
  scale_y_continuous(labels = function(x) sprintf("%.1f Mb", x / 1e6))+
  scale_color_manual(values=md$Color, breaks =md$Group)+
  theme(legend.position = 'top', 
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        legend.key.size = unit(0.3, "lines"))+
  guides(color = guide_legend(nrow = 2))
nx

ggsave('~/artocarpus_pangenome/02_genome_assembly/20260615_Nx_Plot.pdf',nx,height=2.5,width=3)

```

