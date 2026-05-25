#!/bin/bash
#SBATCH -J enrich
#SBATCH -N 1
#SBATCH -n 1  
#SBATCH --cpus-per-task=2
#SBATCH -p qitinglab
#SBATCH -e enrich_%A_%a.err
#SBATCH -o enrich_%A_%a.out
#SBATCH -A qitinglab-grp
#SBATCH --mem=30G
#SBATCH --array=1
#SBATCH --nodelist=qiting-svr1

source ~/.bashrc
conda activate jupyter
enrich_script_dir=/qitinglab/software/select_SNP_enrich
sampleing_py=/qitinglab/software/select_SNP_enrich/python/sampling.py
enrich_py=/qitinglab/software/select_SNP_enrich/python/enrich_analyisis.py
chr=$SLURM_ARRAY_TASK_ID
query_all_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_all
query_lead_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_top
sQTL_cate=('RI_ratio' 'MX_ratio' 'SE_ratio' 'A3_ratio' 'A5_ratio' 'AF_ratio' 'AL_ratio')
sampling_res_dir='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/sampling_res'
### sampling ############
for qtl_type in ${sQTL_cate[@]}
do
QTL_summary=${query_all_dir}/Metabrain_Augment_${qtl_type}_chr${chr}.txt
lead_QTL=${query_lead_dir}/Metabrain_Augment_${qtl_type}_chr${chr}.cis.summary.txt
python ${sampleing_py} --QTL_summary ${QTL_summary} --lead_QTL ${lead_QTL} --output ${sampling_res_dir}/${qtl_type}_chr${chr}_sampling.txt
done
### enrichment ##########
python ${enrich_py} --SNP_annot /qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/per_variant_annotation/all_annot_per_SNP.txt --enrich_list /qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/per_variant_annotation/all_annot_type.txt \
--sampling /qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/sampling_res/merge_sQTL_sampling.txt --SNP_type sQTL_merge --output /qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/enrichment_result/all_sQTL_enrichment.tsv

