# 06_reference_based_ancestry/

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

## Merge outputs:

```bash
#!/bin/bash

#SBATCH --time=1-00:00:00    
#SBATCH --cpus-per-task=48
#SBATCH --mem=256Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

module load miniconda
source activate snp_array

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

if [ ! -f "$OUTDIR/joint.callable.raw.vcf.gz" ]; then
    echo "Filtering vcf"

    # basic filtering
    bcftools view \
        -i 'FILTER="PASS" && F_MISSING < 0.2 && MAF > 0.01' \
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
        plink2 --threads ${t} --vcf "$OUTDIR/joint.filtered.callable.raw.vcf.gz" --keep ${WD}/${subset}.list --allow-extra-chr --set-missing-var-ids @:# \
                --rm-dup --indep-pairwise 50 5 0.1 --maf 0.05 --hwe 1e-10 --max-alleles 2 --min-alleles 2 --out ${WD}/08_plink/${subset}

        #extract, also a vcf
        plink2 --threads ${t} --vcf "$OUTDIR/joint.filtered.callable.raw.vcf.gz" --keep ${WD}/${subset}.list --allow-extra-chr --set-missing-var-ids @:# \
                --extract ${WD}/08_plink/${subset}.prune.in \
                --make-bed --recode vcf bgz --pca --out ${WD}/08_plink/${subset}
        bcftools index --threads ${t} ${WD}/08_plink/${subset}.vcf.gz
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

# create tree
if [ ! -f ${WD}/10_tree/joint.filtered.callable.raw.min4.phy ]; then
    echo "Create tree on snps"
    mkdir -p ${WD}/10_tree
    python ~/symlinks/software/vcf2phylip/vcf2phylip.py -i "$OUTDIR/joint.filtered.callable.raw.vcf.gz" -f --output-folder ${WD}/10_tree
    iqtree --redo -keep-ident -T ${t} -s ${WD}/10_tree/joint.filtered.callable.raw.min4.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000
    iqtree --redo -keep-ident -T ${t} -s ${WD}/10_tree/joint.filtered.callable.raw.min4.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000
fi


```

## DSuite

```bash
#!/bin/bash

#SBATCH --time=1-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate dsuite

TREE=/project/coffea_pangenome/Artocarpus/Comparative_Paper/cafe5/tree.nwk

#run in parallel
/project/coffea_pangenome/Software/Merondun/Dsuite/Build/Dsuite Dtrios --cores 20 -n comp -t ${TREE} dsuite/comp.pop vcfs/Chr01.MM1.vcf.gz vcfs/Chr02.MM1.vcf.gz vcfs/Chr03.MM1.vcf.gz vcfs/Chr04.MM1.vcf.gz vcfs/Chr05.MM1.vcf.gz vcfs/Chr06.MM1.vcf.gz vcfs/Chr07.MM1.vcf.gz vcfs/Chr08.MM1.vcf.gz vcfs/Chr09.MM1.vcf.gz vcfs/Chr10.MM1.vcf.gz vcfs/Chr11.MM1.vcf.gz vcfs/Chr12.MM1.vcf.gz vcfs/Chr13.MM1.vcf.gz vcfs/Chr14.MM1.vcf.gz vcfs/Chr15.MM1.vcf.gz vcfs/Chr16.MM1.vcf.gz vcfs/Chr17.MM1.vcf.gz vcfs/Chr18.MM1.vcf.gz vcfs/Chr19.MM1.vcf.gz vcfs/Chr20.MM1.vcf.gz vcfs/Chr21.MM1.vcf.gz vcfs/Chr22.MM1.vcf.gz vcfs/Chr23.MM1.vcf.gz vcfs/Chr24.MM1.vcf.gz vcfs/Chr25.MM1.vcf.gz vcfs/Chr26.MM1.vcf.gz vcfs/Chr27.MM1.vcf.gz vcfs/Chr28.MM1.vcf.gz

#calculate f4s
/project/coffea_pangenome/Software/Merondun/Dsuite/Build/Dsuite Fbranch ${TREE} dsuite/DTparallel_comp_comp_combined_tree.txt > dsuite/DTparallel_comp_comp_combined_tree.fb

#plot
python3 /project/coffea_pangenome/Software/Merondun/Dsuite/utils/dtools.py dsuite/DTparallel_comp_comp_combined_tree.fb ${TREE} -n dsuite/DSUITE --dpi 300 --ladderize --tree-label-size 10

#get Z scores
/project/coffea_pangenome/Software/Merondun/Dsuite/Build/Dsuite Fbranch ${TREE} dsuite/DTparallel_comp_comp_combined_tree.txt -Z > dsuite/DTparallel_comp_comp_combined_tree.Z
#tail -n 14 dsuite/DTparallel_SamplesN12_N12_combined_tree.Z > dsuite/DTparallel_SamplesN12_N12_combined_tree.Z.txt
 
```

