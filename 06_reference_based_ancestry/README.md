# 06_reference_based_ancestry/

Reference-based approach calling SNPs against HART050 to infer DSuite allele sharing, ADMIXTURE plots, PCA, and a SNP-based phylogeny.

Output:

![ancestry](/imgs/ancestry.png)



___



These scripts use a reference based approach to sample HART050 to call SNPs (including just A. camansi and A. mariannensis as shallow outgroup and A. odoratissimus as a proper root). 

* Subsets HiFi reads, aligns to the reference, calls variants with DeepVariant, and computes depth-based callable regions for each sample.
* Calculates heterozygous SNP counts and per‑base heterozygosity within callable intervals, outputting all summary statistics per sample.

Map hifi, call SNPs:

```bash
#!/bin/bash

#SBATCH --time=4-00:00:00    
#SBATCH --cpus-per-task=6
#SBATCH --mem=48Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

module load apptainer
module load miniconda
#mamba create -n snp_array python bbtools minimap2 samtools bcftools glnexus seqkit r-tidyverse r-ape r-ggtree r-treeio bedtools iqtree admixture bioconductor-snprelate r-ranger r-randomforest r-tidymodels mosdepth 
source activate snp_array

SAMPLE=${1:?Missing SAMPLE argument}
WD=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/reference_based
HIFI=/project/coffea_pangenome/Artocarpus/Concatenated_Reads/${SAMPLE}.HiFi.fastq.gz
GENOME=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/assemblies/unmasked/HART050.fa

t=24

mkdir -p \
  ${WD}/01_subset_reads \
  ${WD}/02_bams \
  ${WD}/03_vcfs \
  ${WD}/04_mosdepth \
  ${WD}/05_callable \
  ${WD}/06_stats

# Subset reads so that we have roughly equal inputs
if [ ! -f ${WD}/01_subset_reads/${SAMPLE}.25gb.fastq.gz ]; then
    echo "Subsetting Reads for: ${SAMPLE}"
    bbduk.sh in=${HIFI} out=${WD}/01_subset_reads/${SAMPLE}.25gb.fastq.gz maxbasesout=25000000000
fi

# Align
if [ ! -f ${WD}/02_bams/${SAMPLE}.sorted.bam ]; then
    echo "Aligning Reads for: ${SAMPLE}"
    minimap2 -ax map-hifi -t ${t} \
        -R @RG\\tID:${SAMPLE}\\tPL:PACBIO\\tLB:${SAMPLE}\\tSM:${SAMPLE} ${GENOME} \
        ${WD}/01_subset_reads/${SAMPLE}.25gb.fastq.gz | \
        samtools view -F 4 -bS - | \
        samtools sort -@ ${t} -o ${WD}/02_bams/${SAMPLE}.sorted.bam
    samtools index ${WD}/02_bams/${SAMPLE}.sorted.bam
fi

# Call variants, just use diploid deep variant model since we will ignore dosage 
# apptainer pull deepvariant_latest.sif docker://google/deepvariant:latest
if [ ! -f ${WD}/03_vcfs/${SAMPLE}.pt.vcf.gz ]; then
    echo "Calling SNPs for: ${SAMPLE}"
    DEEPVAR="/project/coffea_pangenome/Breadfruit_SNP_Array/containers/deepvariant_gh1060.sif"
    apptainer exec \
        -B /project/coffea_pangenome:/project/coffea_pangenome \
        ${DEEPVAR} run_deepvariant \
        --make_examples_extra_args='small_model_call_multiallelics=false' \
        --model_type PACBIO \
        --ref ${GENOME} \
        --reads ${WD}/02_bams/${SAMPLE}.sorted.bam \
        --output_vcf ${WD}/03_vcfs/${SAMPLE}.pt.vcf.gz \
        --output_gvcf ${WD}/03_vcfs/${SAMPLE}.pt.gvcf.gz \
        --sample_name ${SAMPLE} \
        --num_shards ${t} \
        --postprocess_cpus ${t}
    tabix -f -p vcf ${WD}/03_vcfs/${SAMPLE}.pt.vcf.gz
fi

# identify callable regions MQ >=20 dp 0.5x-2x expected 
if [ ! -f ${WD}/04_mosdepth/${SAMPLE}.mosdepth.summary.txt ]; then
    echo "Computing per-base depth (mosdepth) for: ${SAMPLE}"
    mosdepth -t ${t} -Q 20 ${WD}/04_mosdepth/${SAMPLE} ${WD}/02_bams/${SAMPLE}.sorted.bam
fi

# estimate expected depth from mosdepth summary (total/mean)
if [ ! -f ${WD}/06_stats/${SAMPLE}.depth_thresholds.txt ]; then
    echo "Estimating expected depth for ${SAMPLE}"
    MEAN_DP=$(awk '$1=="total"{print $4}' ${WD}/04_mosdepth/${SAMPLE}.mosdepth.summary.txt)

    if [ -z "${MEAN_DP}" ]; then
        echo "ERROR: could not parse mean depth from mosdepth summary: ${WD}/04_mosdepth/${SAMPLE}.mosdepth.summary.txt" >&2
        exit 1
    fi

MIN_DP=$(python3 - <<PY
m=float("${MEAN_DP}")
print(max(1, int(m*0.5)))
PY
)

MAX_DP=$(python3 - <<PY
m=float("${MEAN_DP}")
print(int(m*2.0))
PY
)

    echo "Mean depth=${MEAN_DP}; callable depth range=[${MIN_DP},${MAX_DP}] (MAPQ>=20)" \
    | tee ${WD}/06_stats/${SAMPLE}.depth_thresholds.txt

fi 

# Build callable intervals from per-base depths:
if [ ! -f ${WD}/06_stats/${SAMPLE}.callable_bp.tsv ]; then
    echo "Building callable intervals for ${SAMPLE}"
    CALLABLE_BED=${WD}/05_callable/${SAMPLE}.callable.MQ20.DP${MIN_DP}-${MAX_DP}.bed
    zcat ${WD}/04_mosdepth/${SAMPLE}.per-base.bed.gz | \
    awk -v min=${MIN_DP} -v max=${MAX_DP} 'BEGIN{OFS="\t"} $4>=min && $4<=max {print $1,$2,$3}' | \
    bedtools merge -i - > ${CALLABLE_BED}

    CALLABLE_BP=$(awk '{s+=$3-$2} END{print s+0}' ${CALLABLE_BED})
    echo -e "${SAMPLE}\t${CALLABLE_BP}" > ${WD}/06_stats/${SAMPLE}.callable_bp.tsv
fi 

### heterozygosity time 
if [ ! -f ${WD}/06_stats/${SAMPLE}.output.tsv ]; then
    echo "Counting heterozygous SNPs in callable regions for: ${SAMPLE}"
    VCF=${WD}/03_vcfs/${SAMPLE}.pt.vcf.gz
    MIN_GQ=20
    MIN_VCF_DP=10

    HET_SNPS=$(bcftools view \
    -R ${CALLABLE_BED} \
    -f PASS \
    -v snps \
    -i "GT='het' && GQ>=${MIN_GQ} && DP>=${MIN_VCF_DP}" \
    ${VCF} | \
    bcftools view -H | wc -l)

# heterozygosity per bp (pi-like SNP density proxy)
HET_PER_BP=$(python3 - <<PY
het=int("${HET_SNPS}")
bp=int("${CALLABLE_BP}")
print("nan" if bp==0 else het/bp)
PY
)

    printf "%s\t%s\t%s\t%s\n" \
    "${SAMPLE}" \
    "${CALLABLE_BP}" \
    "${HET_SNPS}" \
    "${HET_PER_BP}" \
    > ${WD}/06_stats/${SAMPLE}.het_density.tsv

    echo "Done: ${SAMPLE} callable_bp=${CALLABLE_BP} het_snps=${HET_SNPS} het_per_bp=${HET_PER_BP}"
    echo -e "${SAMPLE}\t${CALLABLE_BP}\t${HET_SNPS}\t${HET_PER_BP}" > ${WD}/06_stats/${SAMPLE}.output.tsv
fi 
```

## Merge VCFs & Run PopGen

* Joint-call GVCFs with GLnexus, convert to VCF, index, and filter variants to callable regions shared across all samples
* Generate LD‑pruned SNP subsets, convert to PLINK, run ADMIXTURE (K=2–5), and perform PCA.
* Build phylogenetic tree from SNP data (IQ-TREE).

```bash
#!/bin/bash

#SBATCH --time=1-00:00:00    
#SBATCH --cpus-per-task=48
#SBATCH --mem=512Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

#module load miniconda
#source activate snp_array

set -euo pipefail

### check sample names first!
WD=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/reference_based
t=48
REF_SAMPLE="HART050"  
GVCF_DIR=${WD}/03_vcfs
OUTDIR=${WD}/07_jointvcf

mkdir -p ${OUTDIR}

if [ ! -f ${OUTDIR}/joint.raw.vcf.gz ]; then
    echo "Building glnexus db"
    ls ${GVCF_DIR}/*.gvcf.gz > "$OUTDIR/gvcfs.list"

    glnexus_cli \
        --threads 48 \
        --config DeepVariantWGS \
        $(cat "$OUTDIR/gvcfs.list") \
        > "$OUTDIR/joint.bcf"

    bcftools view -Oz \
        -o "$OUTDIR/joint.raw.vcf.gz" \
        "$OUTDIR/joint.bcf"

    bcftools index -t "$OUTDIR/joint.raw.vcf.gz"
fi 

if [ ! -f ${OUTDIR}/joint.filtered.callable.raw.vcf.gz ]; then
    echo "Filtering vcf"

    # basic filtering
    bcftools view \
        -i 'F_MISSING = 0' \
        -Oz \
        -o "$OUTDIR/joint.filtered.vcf.gz" \
        "$OUTDIR/joint.raw.vcf.gz"

    bcftools index -t "$OUTDIR/joint.filtered.vcf.gz"

    # regions callable in all N samples
    # if you want to check 
    for f in ${WD}/05_callable/*.bed; do
        bp=$(awk '{sum += $3 - $2} END {print sum}' "$f")
        prefix=$(basename "$f" .bed)
        echo -e "${prefix}\t${bp}"
    done

    # now filter ensuring all samples overlap
    N=17
    bedtools multiinter \
        -i ${WD}/05_callable/*.callable.MQ20*.bed \
        | awk -v n="$N" '$4 == n {print $1"\t"$2"\t"$3}' \
        | bedtools merge \
        > ${WD}/05_callable/cohort.callable.all${N}.bed

    bp=$(awk '{sum += $3 - $2} END {print sum/1e6}' ${WD}/05_callable/cohort.callable.all${N}.bed)
    echo "There are ${bp} Mb callable bases"

    # Filtered for callable
    bcftools view \
        -R ${WD}/05_callable/cohort.callable.all${N}.bed \
        -Oz \
        -o "$OUTDIR/joint.filtered.callable.raw.vcf.gz" \
        "$OUTDIR/joint.filtered.vcf.gz"

    bcftools index -t "$OUTDIR/joint.filtered.callable.raw.vcf.gz"
fi 

#LD prune, both including camansi / mariannensis and without
mkdir -p ${WD}/08_plink

for subset in Shallow Breadfruit; do 
    if [ ! -f ${WD}/08_plink/${subset}.bim ]; then
        # Extract LD 
        echo "Creating LD pruned subset for ${subset}"
        bcftools view -S ${WD}/${subset}.list ${OUTDIR}/joint.filtered.callable.raw.vcf.gz -m2 -M2 -i 'MAF>0.05' -Ou \
            | bcftools +prune -m 0.1 --window 5kb -Oz -o ${WD}/08_plink/${subset}.pruned.vcf.gz

        # Convert to bed 
        plink2 --threads ${t} --vcf ${WD}/08_plink/${subset}.pruned.vcf.gz \
            -chr-set 28 --allow-extra-chr --set-missing-var-ids @:# \
            --make-bed --recode vcf bgz --out ${WD}/08_plink/${subset}
            
        sed -i 's/^Chr0*\([0-9]\+\)\t/\1\t/' ${WD}/08_plink/${subset}.bim
    fi
done

#run admixture for both subsets , from k 2 to 5 
mkdir -p ${WD}/09_admixture
for subset in Shallow Breadfruit; do 
    for K in {2..5}; do
        if [ ! -f ${WD}/09_admixture/${subset}.log${K}.out ]; then
            cd ${WD}/09_admixture
            echo "Runnign admixture at K: ${K} for ${subset}"
            admixture -j7 --cv=5 ${WD}/08_plink/${subset}.bed ${K} > ${subset}.log${K}.out
        fi
    done
done 

# PCA
if [ ! -f ${WD}/09_admixture/pca_scores.txt ]; then
    echo "Running merothon PCA"
    vcf_to_pca --vcf 08_plink/Shallow.pruned.vcf.gz --out 09_admixture/pca.png
fi

# create tree
if [ ! -f ${WD}/10_tree/joint.filtered.callable.raw.min4.phy.varsites.phy.contree ]; then
    echo "Create tree on snps"
    mkdir -p ${WD}/10_tree
    python ~/symlinks/software/vcf2phylip/vcf2phylip.py -i "$OUTDIR/joint.filtered.callable.raw.vcf.gz" -f --output-folder ${WD}/10_tree
    iqtree --redo -keep-ident -T ${t} -s ${WD}/10_tree/joint.filtered.callable.raw.min4.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000
    iqtree --redo -keep-ident -T ${t} -s ${WD}/10_tree/joint.filtered.callable.raw.min4.phy.varsites.phy --seqtype DNA -m "GTR+ASC" -alrt 1000 -B 1000
fi

```

## Plot ADMIXTURE

* Loads metadata and ADMIXTURE Q‑files, merges them, plots CV‑error and faceted ADMIXTURE barplots for K=2–5.
* Generates labeled PCA plots.

CV Error: combine first:

```
for g in Shallow Breadfruit; do 
	grep "CV" ${g}*out | sed 's/.*log//g' | sed 's/.out:CV//g' > ${g}-ADMIXTURE.CVs.txt
done 
```

Plot:

```R
### Plot ADMIXTURE across landscape (tesselation)
setwd('/project/coffea_pangenome/Artocarpus/Pangenome_Paper/reference_based/09_admixture')
library(tidyverse)
library(meRo) #devtools::install_github('merondun/meRo')
library(ggrepel)

# Hash out which species to run, entire script will run afterwards 
admix_run = 'Breadfruit'
admix_run = 'Shallow'

qdir = '.' #directory with Q files

admix = melt_admixture(prefix = admix_run, qdir = qdir)

#read in metadata
md = read_tsv('~/artocarpus_pangenome/samples.info') %>% dplyr::select(ID = Sample, Name, Group, Color) %>% 
  mutate(Shape = 21) %>% 
  rbind(.,data.frame(ID = c('HART063','HART067','HART061'),
                     Group = c('camansi','mariannensis','odoratissimus'),
                     Name = c('camansi','mariannensis','odoratissimus'),
                     Color = c('black','black','black'),
                     Shape = c(24,25,8)))
md <- md %>% mutate(lab = paste0(Name,'\n(',ID,')'))

admixmd = left_join(admix,md)

# Reorder individuals based on ploidy group 
admixmd = admixmd %>% mutate(ID = factor(ID, levels=unique(ID)),
                             lab = factor(lab, levels=unique(lab)))

# Plot CV error 
cv = read.table(paste0(admix_run,'-ADMIXTURE.CVs.txt'),header=FALSE)
names(cv) <- c('K','d1','d2','Error')
cvs = 
  cv %>% ggplot(aes(x=K,y=Error))+
  geom_line(show.legend = F,col='black')+
  geom_point(show.legend = F,col='black',size=2)+
  scale_color_manual(values=viridis(3))+
  scale_x_continuous(breaks=function(x) pretty(x,n=10))+
  ylab('C-V Error')+
  theme_classic()
cvs
ggsave(paste0('20260615_ADMIXTURE-CV-Error_',admix_run,'.pdf'),cvs,height=2,width=3,dpi=600)

#if you want to add CV error directly on the label
names(cv) = c('Specified_K','d1','d2','Error')
cv = cv %>% mutate(label = paste0('K',Specified_K,' (',round(Error,2),')')) %>% select(!c(d1,d2,Error)) %>% arrange(Specified_K)
cv

# K = 2-5, by haplogroup  
adplot =
  admixmd %>% 
  #filter(Specified_K == 5 | Specified_K == 2) %>%  #specify the levels you want 
  mutate(Specified_K = paste0('K',Specified_K)) %>% 
  ggplot(aes(x = factor(lab), y = Q, fill = factor(K))) +
  geom_col() +
  facet_grid(Specified_K~Group, scales = "free", space = "free") +
  theme_minimal(base_size=6) + labs(x = "",y = "") +
  scale_y_continuous(expand = c(0, 0),n.breaks = 3) +
  scale_fill_manual(values=viridis(5))+
  theme(
    panel.spacing.x = unit(0.01, "lines"),
    #axis.text.x = element_blank(),
    axis.text.x=element_text(angle=90,size=5),
    axis.text.y = element_text(size=3),
    panel.grid = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(size=6)
  )
adplot

ggsave(paste0('20260615_Admixture-K2-5_',admix_run,'.pdf'),adplot,height=2.5,width=2.5,dpi=600)

# Plot PCA
pc <- read_tsv('pca_scores.txt')
eigs <- read_tsv('pca_values.txt')
eigs = eigs %>% mutate(VE = Explained_Variance/sum(Explained_Variance))
pcm <- left_join(pc %>% dplyr::rename(ID = SampleID),md)

# Extract percent variance explained
pc1_pve <- eigs %>% filter(PC == "PC1") %>% pull(Explained_Variance) * 100
pc2_pve <- eigs %>% filter(PC == "PC2") %>% pull(Explained_Variance) * 100

# Plot PCA
pp <- ggplot(pcm, aes(x = PC1, y = PC2, fill = Group, shape = Group)) +
  geom_point(size = 4, stroke = 0.5) +
  # repel labels
  geom_text_repel(
    aes(label = lab),
    size = 2,
    max.overlaps = Inf,
    point.padding = 0.5,
    min.segment.length = 0
  ) +
  
  scale_fill_manual(values = pcm$Color %>% setNames(pcm$Group)) +
  scale_shape_manual(values = pcm$Shape %>% setNames(pcm$Group)) +
  labs(
    x = paste0("PC1 (", round(pc1_pve, 1), "%)"),
    y = paste0("PC2 (", round(pc2_pve, 1), "%)"),
    color = "Group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )
pp

ggsave(paste0('20260615_PCA_',admix_run,'.pdf'),pp,height=5,width=7,dpi=600)


```

## Plot Tree

* Reads and roots the IQ‑TREE phylogeny, plots a labeled SNP tree, saves both fig and a binary ultrametric tree for downstream Dsuite work.

```R
### Plot SNP tree from ref based approach
setwd('/project/coffea_pangenome/Artocarpus/Pangenome_Paper/reference_based/')
library(ggtree)
library(treeio)
library(tidyverse)
library(ape)

md = read_tsv('~/artocarpus_pangenome/samples.info') %>% dplyr::select(ID = Sample, Name, Group, Color) %>% 
  mutate(Shape = 21) %>% 
  rbind(.,data.frame(ID = c('HART063','HART067','HART061'),
                     Group = c('camansi','mariannensis','odoratissimus'),
                     Name = c('camansi','mariannensis','odoratissimus'),
                     Color = c('black','black','black'),
                     Shape = c(24,25,8)))

tree <- '10_tree/joint.filtered.callable.raw.min4.phy.varsites.phy.contree'

iqtree = read.iqtree(tree)
iqtr = root(as.phylo(iqtree),'HART061')
iqtr2 <- drop.tip(iqtr, "HART061")

#plot with outgroups 
md <- md %>% mutate(lab = paste0(Name,' (',ID,')'))
gt <- ggtree(iqtr2, layout = "rectangular") %<+% md
gtp <- gt + 
  geom_nodepoint(mapping=aes(subset=(as.numeric(label) >= 95)),col='black',fill='grey90',pch=23,size=0.75,show.legend=F)+
  geom_tippoint(aes(fill = Group, shape = Group),size=2)+
  geom_tiplab(aes(label = lab),align = FALSE,hjust=-0.25,size=3)+
  scale_fill_manual(values = pcm$Color %>% setNames(pcm$Group)) +
  scale_shape_manual(values = pcm$Shape %>% setNames(pcm$Group)) +  guides(fill=guide_legend(nrow=4,override.aes=list(shape=21)),
                                                                           shape=guide_legend(nrow=5))+
  theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10),
        legend.key.size = unit(0.2, "cm"),    legend.position = 'top') + 
  xlim(c(min(gt$data$x,na.rm=TRUE),max(gt$data$x,na.rm=TRUE)*1.3))
gtp

ggsave(paste0('20260615_Tree.pdf'),gtp,height=4,width=3.5,dpi=600)

# Save tree ultrametric with resolved polytomies for dsuite
simp <- iqtr
simp$node.label <- NULL
rutree <- multi2di(simp, random = TRUE)
rutree_d = root(rutree,'HART061',resolve.root = TRUE)
plot(rutree_d, cex = 0.6)
is.rooted(rutree_d)
is.binary(rutree_d)
is.ultrametric(rutree_d)
write.tree(rutree_d,'/project/coffea_pangenome/Artocarpus/Pangenome_Paper/reference_based/binary_tree.nwk')
```



## DSuite

* Extracts per‑chromosome VCFs, runs Dsuite (DtriosParallel, Fbranch, Z‑scores) using the binary rooted SNP tree. 

```bash
#!/bin/bash

#SBATCH --time=0-06:00:00    
#SBATCH --cpus-per-task=48
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

#module load miniconda
#source activate snp_array

set -euo pipefail

### check sample names first!
WD=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/reference_based
t=48

# Dsuite -t ${WD}/10_tree/joint.filtered.callable.raw.min4.phy.varsites
if [ ! -f ${WD}/11_dsuite/DTparallel_SETS_breadfruit_combined_tree.fb ]; then
    echo "Running dsuite"
    mkdir -p ${WD}/11_dsuite/chrs
    TREE=/project/coffea_pangenome/Artocarpus/Pangenome_Paper/reference_based/binary_tree.nwk

    for chr in $(bcftools query -f '%CHROM\n' ${OUTDIR}/joint.filtered.callable.raw.vcf.gz | sort -u); do
        echo "Extracting $chr"
        bcftools view -r ${chr} \
            -Oz -o ${WD}/11_dsuite/chrs/${chr}.vcf.gz \
            ${OUTDIR}/joint.filtered.callable.raw.vcf.gz
        tabix -p vcf ${WD}/11_dsuite/chrs/${chr}.vcf.gz
    done

    # run dsuite
    ~/symlinks/software/Dsuite/utils/DtriosParallel -t ${TREE} -n breadfruit ${WD}/11_dsuite/SETS.txt \
        ${WD}/11_dsuite/chrs/*vcf.gz

    #calculate f4s
    ~/symlinks/software/Dsuite/Build/Dsuite Fbranch ${TREE} ${WD}/11_dsuite/DTparallel_SETS_breadfruit_combined_tree.txt > ${WD}/11_dsuite/DTparallel_SETS_breadfruit_combined_tree.fb

    #plot
    python3 ~/symlinks/software/Dsuite/utils/dtools.py --outgroup HART061 ${WD}/11_dsuite/DTparallel_SETS_breadfruit_combined_tree.fb ${TREE} -n ${WD}/11_dsuite/breadfruit --dpi 300 --ladderize --tree-label-size 10

    #get Z scores
    ~/symlinks/software/Dsuite/Build/Dsuite Fbranch ${TREE} ${WD}/11_dsuite/DTparallel_SETS_breadfruit_combined_tree.txt -Z > ${WD}/11_dsuite/DTparallel_SETS_breadfruit_combined_tree.Z
    #tail -n 14 dsuite/DTparallel_SamplesN12_N12_combined_tree.Z > dsuite/DTparallel_SamplesN12_N12_combined_tree.Z.txt

fi


```

## ADMIXTOOLS

F4 statistics:

```R
#install.packages("devtools") # if "devtools" is not installed already
#devtools::install_github("uqrmaie1/admixtools")
setwd('/project/coffea_pangenome/Artocarpus/Pangenome_Paper/reference_based/08_plink')
library(admixtools)
library(tidyverse)

prefix <- "Full"

fam <- read_table(
  paste0(prefix, ".fam"),
  col_names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"),
  show_col_types = FALSE
)

fam_fixed <- fam %>%
  mutate(FID = IID)

write.table(
  fam_fixed,
  paste0(prefix, ".fam"),col.names = F, row.names = F, quote=F)

cultivars <- c(
  "HART053", "HART046", "HART049", "HART050",
  "ZZ3", "H6", "HART033", "ZZ7", "ZZ9",
  "HART032", "HART038", "HART030", "HART001", "HART069"
)

pops <- c("HART061", "HART063", "HART067", cultivars)

extract_f2(
  prefix,
  my_f2_dir,
  pops = pops,
  auto_only = FALSE,
  maxmiss = 1,
  overwrite = TRUE
)

f2_blocks <- f2_from_precomp(my_f2_dir, pops = pops)

d_res <- list()

for (cultivar in cultivars) {
  cat("Working on:", cultivar, "\n")
  
  d_res[[cultivar]] <- f4(
    f2_blocks,
    pop1 = "HART063",
    pop2 = "HART067",
    pop3 = cultivar,
    pop4 = "HART061",
    f4mode = FALSE
  ) %>%
    mutate(
      cultivar = cultivar,
      ci_low = est - 1.96 * se,
      ci_high = est + 1.96 * se,
      interpretation = case_when(
        est > 0 ~ "mariannensis-like",
        est < 0 ~ "camansi-like",
        TRUE ~ "balanced"
      )
    )
}

d_res <- bind_rows(d_res) %>%
  arrange(est)

write_tsv(d_res, "D_HART063_HART067_cultivar_HART061.tsv")

}
d_res <- map_dfr(
  cultivars,
  \(cultivar) {
    f4(
      prefix,
      pop1 = "HART063",
      pop2 = "HART067",
      pop3 = cultivar,
      pop4 = "HART061",
      f4mode = FALSE
    ) %>%
      mutate(cultivar = cultivar)
  }
) %>%
  mutate(
    ci_low = est - 1.96 * se,
    ci_high = est + 1.96 * se,
    interpretation = case_when(
      est > 0 ~ "mariannensis-like",
      est < 0 ~ "camansi-like",
      TRUE ~ "balanced"
    )
  ) %>%
  arrange(est)

write_tsv(d_res, "D_HART063_HART067_cultivar_HART061.tsv")

ggplot(
  d_res,
  aes(
    x = est,
    y = fct_reorder(cultivar, est)
  )
) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0) +
  geom_point(size = 2.5) +
  labs(
    x = "D(HART063, HART067; cultivar, HART061)",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "D_HART063_HART067_cultivar_HART061.pdf",
  width = 6,
  height = 4.5
)
```

