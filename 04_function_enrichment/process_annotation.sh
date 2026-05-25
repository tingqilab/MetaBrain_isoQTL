#!/bin/sh
#SBATCH -J make_snp_annotation
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH -p qitinglab
#SBATCH -e make_snp_annotation_%A_%a.err
#SBATCH -o make_snp_annotation_%A_%a.out
#SBATCH -A qitinglab-grp
#SBATCH --mem=2G
#SBATCH --nodelist=qiting-svr1
#SBATCH --array=5001-6336%50
## 1-1000 yes
## 1001-2000 yes
## 2001-3000 yse
## 3001-4000
## 4001-5000
## 5001-6000
## 6001-6336


########### annotaion filter #####################################################################

## RBP filter##
# extract significant peaks #############
# RBP_bed_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/function_annotation/RBP
# RBP_clean_bed_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/function_annotation_clean/RBP
# bed_list=(${RBP_bed_dir}/*.bed)
# bed=${bed_list[$SLURM_ARRAY_TASK_ID]}
# basename=$(basename ${bed} _OneFold_1e-05Pvalue.bed)
# awk -F '\t' -v name="$basename" 'BEGIN{OFS="\t"} {print $0, name}' "$bed" > "${RBP_clean_bed_dir}/${basename}.bed"

# Histone filter ##
# histone_bed='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/function_annotation/Histon/RoadMap_sample'
# histone_bed_clean='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/function_annotation_clean/Histon'
# bed_list=(${histone_bed}/*.bed)
# bed=${bed_list[$SLURM_ARRAY_TASK_ID]}
# basename=$(basename ${bed} .bed)
# ## extracting significant peak
# awk '{if($7>1 && $8>5)print}' ${bed} | awk -F '\t' {'printf ("%s\t%s\t%s\t%s\n", $1, $2, $3, $4 )'} >  ${histone_bed_clean}/${basename}_OneFold_1e-05Pvalue.bed

# cCRE filter ##
# cCRE_bed=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/function_annotation/cCRE/cCRE_brain
# cCRE_bed_clean=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/function_annotation_clean/cCRE
# bed_list=(${cCRE_bed}/*.bed)
# bed=${bed_list[$SLURM_ARRAY_TASK_ID]}
# awk -F '\t' {'printf ("%s\t%s\t%s\t%s\n", $1, $2, $3, $10 )'} ${bed} > ${cCRE_bed_clean}/control_brain$((SLURM_ARRAY_TASK_ID+1)).bed


## TF filter, CTCF peak
# TF_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/function_annotation/TF
# TF_clean_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/function_annotation_clean/TF
# bed_gz_list=(${TF_dir}/*.bed.gz)
# bed_gz=${bed_gz_list[$SLURM_ARRAY_TASK_ID]}
# basename=$(basename ${bed_gz} .bed.gz)
# zcat ${bed_gz} | awk -F '\t' {'printf ("%s\t%s\t%s\n", $1, $2, $3)'} > ${TF_clean_dir}/${basename}.bed


############################## intersetion variant and these annotation ####################################################
############################## cCRE, chromatin_stats ,Histon , RBP, TF ###########################
### bedtools intersect -wa -wb #######################
### use this annotation to annotate all of the snp ################

# source ~/.bashrc
# conda activate flair

# # === 输入路径 ===
# SNP_bed_dir='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/all_variant/bed_file'
# annotation_bed_file_list='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/function_annotation_clean/annot_bed.list'
# out_dir='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/variant_function_intersection'

# num_annot=288
# num_chr=22
# task_id=$SLURM_ARRAY_TASK_ID
# annot_idx=$(( (task_id - 1) / num_chr + 1 ))
# chr_idx=$(( (task_id - 1) % num_chr + 1 ))

# # === 获取 annotation 文件路径 ===
# annot_bed=$(sed -n "${annot_idx}p" "$annotation_bed_file_list")
# annot_name=$(basename "$annot_bed" .bed)
# annot_type=$(basename $(dirname "$annot_bed"))
# annotation_base=${annot_type}_${annot_name}

# chr_name="chr${chr_idx}"
# snp_bed=${SNP_bed_dir}/${chr_name}.bed
# mkdir -p ${out_dir}/${chr_name}_annot
# save_dir=${out_dir}/${chr_name}_annot

# # === 执行交集 ===
# echo "[$(date)] intersecting variant ${snp_bed} x ${annotation_base}"
# bedtools intersect -wa -wb \
#   -a ${snp_bed} \
#   -b ${annot_bed} \
#   > ${save_dir}/${annotation_base}.txt