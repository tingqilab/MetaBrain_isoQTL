#!/bin/sh

work_dir=03_metabrain_cohort_data
cohorts=("BM10" "BM22" "BM36" "BM44" "ROSMAP" "BrainGVEx" "CMC" "HBCC" "LIBD" "Mayo")
cohort=${cohorts[$SLURM_ARRAY_TASK_ID-1]}
hg38_genotype_dir=${work_dir}/${cohort}/genotype_hg38
cd ${hg38_genotype_dir}
${plink} --bfile ./${cohort}_hg38_merged --make-bed --out ./${cohort}_hg38_merged_sort
${plink} --bfile ./${cohort}_hg38_merged_sort --indep-pairwise 50 5 0.2 --out ./${cohort}_hg38_merged_pruned_data
${plink} --bfile ./${cohort}_hg38_merged_sort --extract ./${cohort}_hg38_merged_pruned_data.prune.in --make-bed --out ./${cohort}_hg38_genotype_pruned
${plink} --bfile ./${cohort}_hg38_genotype_pruned --pca 10 --out ./${cohort}_hg38_genotype_pca
