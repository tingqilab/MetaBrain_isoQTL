#!/bin/sh
#SBATCH -J make_annot
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH -p qitinglab
#SBATCH -o make_annot_%A_%a.out
#SBATCH -e make_annot_%A_%a.err
#SBATCH -A qitinglab-grp
#SBATCH --mem=5G
#SBATCH --nodelist=qiting-svr2
#SBATCH --array=1-22

# source ~/.bashrc
# conda activate jupyter

# 脚本路径
# bulid_index_py=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_trait_hertability_enrichment/SLDSC/script/bulid_index_for_hm3QTL.py

# # QTL types
# QTL_types=('isoQTL' 'irQTL')
# num_chr=22

# # 计算当前任务的 QTL_type 和 chr
# qtl_idx=$(( (SLURM_ARRAY_TASK_ID - 1) / num_chr ))
# chr=$(( (SLURM_ARRAY_TASK_ID - 1) % num_chr + 1 ))
# QTLtype=${QTL_types[$qtl_idx]}

# echo "Task ID: $SLURM_ARRAY_TASK_ID"
# echo "QTL type: $QTLtype"
# echo "Chromosome: $chr"

# # 调用 Python 脚本
# python ${bulid_index_py} --QTLtype $QTLtype --chr $chr
############### 2.index the sampling pannel ###################
##### hm3 index and general index
##### hm3 index for ldsc
##### general index for function enrichment
# 1) make the index file
# index_sampling_panel_py=/qitinglab/shenlongxin/metabrain/ANALYSIS_PROCESS/4_trait_heritability/SLDSC/script/index_sampling_panel.py
# echo "process eQTL Gencode"
# python ${index_sampling_panel_py} /qitinglab/shenlongxin/metabrain/ANALYSIS_PROCESS/4_trait_heritability/SLDSC/sampling_pannel eQTL Gencode
# echo "process sQTL Gencode"
# python ${index_sampling_panel_py} /qitinglab/shenlongxin/metabrain/ANALYSIS_PROCESS/4_trait_heritability/SLDSC/sampling_pannel sQTL Gencode
# echo "process eQTL Augment"
# python ${index_sampling_panel_py} /qitinglab/shenlongxin/metabrain/ANALYSIS_PROCESS/4_trait_heritability/SLDSC/sampling_pannel eQTL Augment
# echo "process sQTL Augment"
#python ${index_sampling_panel_py} /qitinglab/shenlongxin/metabrain/ANALYSIS_PROCESS/4_trait_heritability/SLDSC/sampling_pannel sQTL Augment
# 2) create readable index

############### 2.index the sampling pannel ###################



################# 3. make ldsc annot file ###########################
source ~/.bashrc
conda activate jupyter
CHR=$SLURM_ARRAY_TASK_ID
echo "Processing chr${CHR}..."
pycode=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_trait_hertability_enrichment/SLDSC/script/make_annot_file.py
python ${pycode} --chr ${CHR} --QTL_type sQTL
python ${pycode} --chr ${CHR} --QTL_type eQTL
python ${pycode} --chr ${CHR} --QTL_type known_irQTL
python ${pycode} --chr ${CHR} --QTL_type novel_irQTL
python ${pycode} --chr ${CHR} --QTL_type irQTL
python ${pycode} --chr ${CHR} --QTL_type known_isoQTL
python ${pycode} --chr ${CHR} --QTL_type novel_isoQTL
python ${pycode} --chr ${CHR} --QTL_type isoQTL

