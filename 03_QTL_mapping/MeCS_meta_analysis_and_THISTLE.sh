#!/bin/sh
#SBATCH -J MeCS
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH -p qitinglab
#SBATCH -e MeCS_meta_%A_%a.err
#SBATCH -o MeCS_meta_%A_%a.out
#SBATCH -A qitinglab-grp
#SBATCH --mem=20G
#SBATCH --array=1-462


# ===  MeCS for meta-analysis ##########
###################### 01. make the meta list for meta analysis ####################################
beqtl_dir_eQTL=03_metabrain_cohort_data/Meta_analysis/cohort_besd_format
MeCS_list=03_metabrain_cohort_data/Meta_analysis/cohort_besd_format/meta_besd_list

# 所有cohort名称
cohorts=("CMC" "BrainGVEx" "Mayo" "BM10" "BM22" "BM36" "BM44" "LIBD" "HBCC" "ROSMAP")
QTL_types=("A3_ratio" "A5_ratio" "AL_ratio" "AF_ratio" "MX_ratio" "RI_ratio" "SE_ratio" "transcript_count" "gene_count" "transcript_ratio")

# 保证输出目录存在
mkdir -p ${MeCS_list}
for QTL_type in ${QTL_types[@]};do
    # ------------------ 主循环 ------------------
    for chr in {1..22}; do

        # 初始化四个list
        Gencode_eQTL_list=${MeCS_list}/Gencode_${QTL_type}_chr${chr}.list
        Metabrain_eQTL_list=${MeCS_list}/Metabrain_${QTL_type}_chr${chr}.list
        Metabrain_Augment_eQTL_list=${MeCS_list}/Metabrain_Augment_${QTL_type}_chr${chr}.list
        > ${Gencode_eQTL_list}
        > ${Metabrain_eQTL_list}
        > ${Metabrain_Augment_eQTL_list}

        
        # 遍历cohort
        for cohort in "${cohorts[@]}"; do

            echo "03_metabrain_cohort_data/Meta_analysis/cohort_besd_format/Gencode/${cohort}_${QTL_type}_chr${chr}" >> ${Gencode_eQTL_list}
            echo "03_metabrain_cohort_data/Meta_analysis/cohort_besd_format/Metabrain/${cohort}_${QTL_type}_chr${chr}" >> ${Metabrain_eQTL_list}
            echo "03_metabrain_cohort_data/Meta_analysis/cohort_besd_format/Metabrain_Augment/${cohort}_${QTL_type}_chr${chr}" >> ${Metabrain_Augment_eQTL_list}
        done
        echo "[chr${chr} ${QTL_type}] list files generated."
    done
done
echo "✅ All eQTL and isoQTL list files generated under: ${MeCS_list}"

############################ 02. perform meta analysis with MeCS #####################################
MeCS_list_file=03_metabrain_cohort_data/Meta_analysis/cohort_besd_format/MeCS_psi.list
MeCS_output_BESD=03_metabrain_cohort_data/Meta_analysis/Meta_QTL/besd
MeCS_output_query_all=03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_all
MeCS_output_query_summary=03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_top
MeCS_output_query_significant=03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_5e8
MeCS_output_query_weak_significant=03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_5e6

LISTS=($(cat ${MeCS_list_file}))
TARGET=${LISTS[$SLURM_ARRAY_TASK_ID-1]}
basename=$(basename ${TARGET} .list)

echo "[$(date)] Start processing: ${basename}"

# Step 1: meta-analysis
if [[ -f ${MeCS_output_BESD}/${basename}.besd ]]; then
    echo "✅ Step 1 skipped: ${basename}.besd already exists."
else
    echo "▶️  Step 1 running: MeCS meta-analysis..."
    ${OSCA} --besd-flist ${TARGET} --mecs --out ${MeCS_output_BESD}/${basename} --thread-num 4
    # Step 2: descriptive summary
    ${SMR} --beqtl-summary ${MeCS_output_BESD}/${basename} --descriptive-cis --out ${MeCS_output_query_summary}/${basename}

    # Step 3: all query
    ${SMR} --beqtl-summary ${MeCS_output_BESD}/${basename} --query 1 --out ${MeCS_output_query_all}/${basename}

    # Step 4: significant query
    ${SMR} --beqtl-summary ${MeCS_output_BESD}/${basename} --query 5e-08 --out ${MeCS_output_query_significant}/${basename}
    
    # Step 5: weak significant query
    ${SMR} --beqtl-summary ${MeCS_output_BESD}/${basename} --query 5e-06 --out ${MeCS_output_query_weak_significant}/${basename}
    echo "[$(date)] Finished ${basename}"
fi

########################### 03.THISTLE analysis Using MeCS Brain-meta isoform-eQTL results ####################################

MeCS_output_BESD=03_metabrain_cohort_data/Meta_analysis/Meta_QTL/besd
MeCS_output_query_all=03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_all
MeCS_output_query_summary=03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_top
MeCS_output_query_significant=03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_5e8
MeCS_output_query_weak_significant=03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_5e6
annotation_bed=03_metabrain_cohort_data/Augment_gene.bed
### get this file
LIST_FILE=03_metabrain_cohort_data/Meta_analysis/Meta_QTL/isoQTL.list
### /qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/Meta_analysis/Meta_QTL/besd/Gencode_isoQTL_chr18.besd
### basename=Gencode_isoQTL_chr18
# 获取当前任务对应的目标文件
LISTS=($(cat ${LIST_FILE}))
TARGET=${LISTS[$SLURM_ARRAY_TASK_ID-1]}
basename=$(basename ${TARGET})
# 将 isoQTL 改成 sQTL
basename_sQTL=$(echo ${basename} | sed 's/\(.*\)_isoQTL_\(.*\)/\1_sQTL_\2/')
echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: ${basename} -> ${basename_sQTL}"
# Step 1: sQTL mapping with THISTLE
# if [ ! -f "${MeCS_output_BESD}/${basename_sQTL}.besd" ]; then
${OSCA} --sqtl --beqtl-summary ${TARGET} \
    --bed ${annotation_bed} \
    --thread-num 10 \
    --task-num 1 \
    --task-id 1 \
    --to-smr \
    --out ${MeCS_output_BESD}/${basename_sQTL}
${SMR} --beqtl-summary ${MeCS_output_BESD}/${basename_sQTL} --descriptive-cis --out ${MeCS_output_query_summary}/${basename_sQTL}
${SMR} --beqtl-summary ${MeCS_output_BESD}/${basename_sQTL} --query 1 --out ${MeCS_output_query_all}/${basename_sQTL}
${SMR} --beqtl-summary ${MeCS_output_BESD}/${basename_sQTL} --query 5e-08 --out ${MeCS_output_query_significant}/${basename_sQTL}
${SMR} --beqtl-summary ${MeCS_output_BESD}/${basename_sQTL} --query 5e-06 --out ${MeCS_output_query_weak_significant}/${basename_sQTL} 
else
    echo "Step 1 skipped: BESD file already exists."
fi
