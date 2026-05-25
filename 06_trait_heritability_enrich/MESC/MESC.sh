#!/bin/sh
#SBATCH -J MESC
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH -p qitinglab
#SBATCH -e MESC_%A_%a.err
#SBATCH -o MESC_%A_%a.out
#SBATCH -A qitinglab-grp
#SBATCH --mem=10G
#SBATCH --array=1-112%30


# 1. maniputating the QTL summary format for MESC needing
# source ~/.bashrc
# conda activate jupyter

# query_all_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_all
# mesc_qtl_format_out_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_trait_hertability_enrichment/MESC/QTL_mesc_format
# PYcode='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_trait_hertability_enrichment/MESC/QTL_format_mesc.py'
# QTL_types=('known_isoQTL' 'novel_isoQTL' 'known_irQTL' 'novel_irQTL')
# chr_list=($(seq 1 22))
# qtl_index=$(( SLURM_ARRAY_TASK_ID / 22 ))
# chr_index=$(( SLURM_ARRAY_TASK_ID % 22))
# QTL_type=${QTL_types[$qtl_index]}
# CHR=${chr_list[$chr_index]}
# prefix=${query_all_dir}/Metabrain_Augment_${QTL_type}
# # === running ===
# echo "Running MESC job for prefix=${prefix}"
# python ${PYcode} --N 2865 --chr ${CHR} --prefix ${prefix}

### 2. MESC estimate the expression score(gene set (seperating novel transcript, known transcript))
# source ~/.bashrc
# conda activate mesc

# MESC_DIRT="/qitinglab/software/mesc"
# BASE_DIRT="/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_trait_hertability_enrichment/MESC"

# # ------------------------------
# # 2️⃣ 参数组合定义
# # ------------------------------
# QTL_types=('known_isoQTL' 'novel_isoQTL' 'known_irQTL' 'novel_irQTL')
# annot_type="Metabrain_Augment"
# total_chr=22
# task_id=${SLURM_ARRAY_TASK_ID}   # 0-65
# qtl_idx=$(( task_id / total_chr ))    # 0,1,2
# chr=$(( task_id % total_chr + 1 ))    # 1-22
# QTL_type=${QTL_types[$qtl_idx]}


# echo "--------------------------------------------"
# echo "Task ID: $task_id"
# echo "Chromosome: $chr"
# echo "QTL type: $QTL_type"
# echo "Annotation type: $annot_type"
# echo "--------------------------------------------"


# python ${MESC_DIRT}/run_mesc.py \
#   --compute-expscore-sumstat \
#   --eqtl-sumstat ${BASE_DIRT}/QTL_mesc_format/${annot_type}_${QTL_type}_chr${chr}.txt \
#   --out ${BASE_DIRT}/expre_score/${annot_type}_${QTL_type}


# 3.Estimation of h2med
# ## 1). converting GWAS to sumstats format frist
# ### --h2med [SUMSTAT_FILE] specifies the GWAS summary statistic file in .sumstats file format. We recommend using the munge_sumstats.py script (part of LD Score Regression) to convert summary statistics to this format.

## 2). calculate h2med #####
source ~/.bashrc
conda activate mesc

MESC_DIRT="/qitinglab/software/mesc"
GWAS_DIRT="/qitinglab/sharedata/GWAS_summary/brain_disease_gwas_newest_EUR/sumstats"
BASE_DIRT="/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_trait_hertability_enrichment/MESC"

# 1-based 索引，空字符串占位
# 定义列表
PHENO=("" "AD_2022_NatGenet" "ADHD2022_NatGenet" "age_at_menarche" \
  "ALS2021" "Brain_volume_adjHeight_2020_NC" "EA5_NatHumBehav_2024_MAMA" \
  "Insomnia" "IQ_NG_2018" "neuroticism_new" "PD_Nalls_FinnGen_META" \
  "PGC3_SCZ_wave3.european.autosome.public.v3" "PGC_ASD_2019_NatGenet" \
  "PGC_BIP2024_eur_no23andMe" "PGC_MDD_2025_no23andMe_eur" \
  "PGC_PTSD_2024NG_eur" "SmokingInitiation" )
TRAIT_NAME=("" "AD" "ADHD" "AAM" "ALS" "BV" "EA" "Insomnia" "IQ" \
  "Neuroticism" "PD" "SCZ" "ASD" "BIP" "MDD" "PTSD" "SmokingInit")
methods=("Metabrain_Augment_eQTL" "Metabrain_Augment_isoQTL" "Metabrain_Augment_irQTL" "Metabrain_Augment_known_irQTL" "Metabrain_Augment_known_isoQTL" "Metabrain_Augment_novel_irQTL" "Metabrain_Augment_novel_isoQTL")
# -----------------------------
# 计算 method 和 trait 的索引
# -----------------------------
TASK_ID=$SLURM_ARRAY_TASK_ID

num_traits=16
method_idx=$(( (TASK_ID - 1) / num_traits ))
trait_idx=$(( (TASK_ID - 1) % num_traits + 1 ))

method=${methods[$method_idx]}
j=$trait_idx

echo "--------------------------------------"
echo "Task ID: $TASK_ID"
echo "Method: $method"
echo "Trait: ${TRAIT_NAME[$j]}"
echo "--------------------------------------"

# -----------------------------
# 执行 h2med 计算
# -----------------------------
python ${MESC_DIRT}/run_mesc.py \
    --h2med ${GWAS_DIRT}/${PHENO[$j]}.sumstats.gz \
    --exp-chr ${BASE_DIRT}/expre_score/${method} \
    --out ${BASE_DIRT}/result/${method}_${TRAIT_NAME[$j]}





### LDSC ####
# sumstats_out='/qitinglab/sharedata/GWAS_summary/brain_related/sumstats'
# conda activate LDSC
# ldsc_sumstats='/home/chenhuilin/software/ldsc_software/munge_sumstats.py'
# GWAS_dir=/qitinglab/sharedata/GWAS_summary/brain_related
# for file in ${GWAS_dir}/*.txt
# do
# basename=$(basename $file .txt)
# python ${ldsc_sumstats} --sumstats ${file} --out ${sumstats_out}/${basename}
# done
### LDSC compelet###

## LDSC ####
# sumstats_out='/qitinglab/sharedata/GWAS_summary/brain_disease_gwas_newest_EUR/sumstats'
# ldsc_sumstats='/qitinglab/software/ldsc_software/munge_sumstats.py'
# GWAS=/qitinglab/sharedata/GWAS_summary/brain_disease_gwas_newest_EUR/EA5_NatHumBehav_2024_MAMA.txt
# basename=$(basename $GWAS .txt)
# python ${ldsc_sumstats} --sumstats ${GWAS} --out ${sumstats_out}/${basename}
# ## LDSC compelet###