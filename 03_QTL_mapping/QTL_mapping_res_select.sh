#!/bin/bash

## 1. TPM based phenotype filtering
cohorts=("CMC" "BM10" "BM22" "BM36" "BM44" "ROSMAP" "BrainGVEx" "HBCC" "LIBD" "Mayo")
cohort=${cohorts[$SLURM_ARRAY_TASK_ID-1]}
py_script=03_metabrain_cohort_data/script/TPM_based_phenotype_filter.py
python ${py_script} --cohort ${cohort}

## 2. make bed file
cohorts=(BM10 BM22 BM36 BM44 HBCC LIBD Mayo CMC ROSMAP BrainGVEx)
annotations=(Gencode Metabrain_Augment)
feature_types=(gene transcript)
measure_types=(abundance_count correct_abundance correct_count)

N_COHORT=${#cohorts[@]}        # 10
N_ANNOT=${#annotations[@]}     # 3
N_FEATURE=${#feature_types[@]} # 2
N_MEASURE=${#measure_types[@]} # 4

task_id=${SLURM_ARRAY_TASK_ID}

measure_idx=$(( task_id % N_MEASURE ))
feature_idx=$(( (task_id / N_MEASURE) % N_FEATURE ))
annot_idx=$(( (task_id / (N_MEASURE * N_FEATURE)) % N_ANNOT ))
cohort_idx=$(( task_id / (N_MEASURE * N_FEATURE * N_ANNOT) ))

cohort=${cohorts[$cohort_idx]}
annotation=${annotations[$annot_idx]}
feature=${feature_types[$feature_idx]}
measure=${measure_types[$measure_idx]}

echo "Task ID      : ${task_id}"
echo "Cohort       : ${cohort}"
echo "Annotation   : ${annotation}"
echo "Feature type : ${feature}"
echo "Measure type : ${measure}"

py_script=03_metabrain_cohort_data/script/make_phenotype_bed_file.py

python ${py_script} \
  --cohort ${cohort} \
  --annotation ${annotation} \
  --feature ${feature} \
  --measure ${measure}

## 3. calculate PCA-based hidden covariates and regressed PC and known cov to get residuals for QTL mapping
cohorts=("HBCC" "LIBD" "BM10" "BM22" "BM36" "BM44" "ROSMAP" "BrainGVEx" "CMC" "Mayo")
annotations=("Gencode" "Metabrain" "Metabrain_Augment")
cohort_idx=$(( SLURM_ARRAY_TASK_ID % 10 ))
annotation_idx=$(( SLURM_ARRAY_TASK_ID / 10))
cohort=${cohorts[$cohort_idx]}
annotation=${annotations[$annotation_idx]}
R_script=03_metabrain_cohort_data/script/PCA_based_phenotype_process.R
Rscript ${R_script} --cohort ${cohort} --annotation ${annotation} --tier=-20,-15,-10,-5,0,5,10 --threads 4

### 4.split bed to different chr
cohorts=("HBCC" "LIBD" "BM10" "BM22" "BM36" "BM44" "ROSMAP" "BrainGVEx" "CMC" "Mayo")
annotations=("Gencode" "Metabrain" "Metabrain_Augment")
cohort_idx=$(( SLURM_ARRAY_TASK_ID % 1 ))
annotation_idx=$(( SLURM_ARRAY_TASK_ID / 1))
cohort=${cohorts[$cohort_idx]}
annotation=${annotations[$annotation_idx]}
Residual_dir=03_metabrain_cohort_data/ROSMAP/splicing_phenotype/phenotype_Residual
split_script=03_metabrain_cohort_data/script/split_bedfile.py
python ${split_script} --Residual ${Residual_dir}

# 5. make isoform ratio phenotype bed
make_isoform_ratio_phenotype_script=/03_metabrain_cohort_data/script/make_isoform_ratio_from_transcript_abundance.py
python ${make_isoform_ratio_phenotype_script} --cohort ${cohort} --annotation ${annotation} 


# PCA for splicing phenotype, such as local splice events and isoform ratio
work_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data
PCA_R_script=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/script/PCA_based_phenotype_process_normal_matrix.R
cohorts=("HBCC" "LIBD" "BM10" "BM22" "BM36" "BM44" "ROSMAP" "BrainGVEx" "CMC" "Mayo")
annotations=("Gencode" "Metabrain_Augment")
splicing_phs=("A3_ratio_new" "A5_ratio_new" "AF_ratio_new" "AL_ratio_new" "MX_ratio_new" "RI_ratio_new" "SE_ratio_new")
num_cohorts=1
num_annotations=3
num_splicing=7
cohort_idx=$(( SLURM_ARRAY_TASK_ID % num_cohorts ))
annotation_idx=$(( (SLURM_ARRAY_TASK_ID / num_cohorts) % num_annotations ))
splicing_idx=$(( SLURM_ARRAY_TASK_ID / (num_cohorts * num_annotations) ))
cohort=${cohorts[$cohort_idx]}
annotation=${annotations[$annotation_idx]}
splicing_ph=${splicing_phs[$splicing_idx]}
splicing_ph_prefix=${splicing_ph/_new/}
bed_file=${work_dir}/${cohort}/splicing_phenotype/filtered_result_bed/${cohort}_${annotation}_${splicing_ph}.bed
cov_file=${work_dir}/${cohort}/covariates/${cohort}_known_biological_covariates.txt
mkdir -p ${work_dir}/${cohort}/splicing_phenotype/phenotype_Residual
output_dir=${work_dir}/${cohort}/splicing_phenotype/phenotype_Residual
prefix=${cohort}_${annotation}_${splicing_ph_prefix}

Rscript ${PCA_R_script} \
  --bed ${bed_file} \
  --cov ${cov_file} \
  --tier=-10,-5,0,5,10,15,20 --threads 3 \
  --normalization 'ratio_norma' \
  --output ${output_dir} --prefix ${prefix}


# get the best PC and optimal QTL result for multiply tier output
source ~/.bashrc
conda activate jupyter
py_script=03_metabrain_cohort_data/script/select_save_best_PCtier_res.py
## parameter
cohorts=("HBCC" "LIBD" "BM10" "BM22" "BM36" "BM44" "ROSMAP" "BrainGVEx" "CMC" "Mayo")
annotations=("Gencode" "Metabrain" "Metabrain_Augment")
cohort_idx=$(( SLURM_ARRAY_TASK_ID % 10 ))
annotation_idx=$(( SLURM_ARRAY_TASK_ID / 10))
cohort=${cohorts[$cohort_idx]}
annotation=${annotations[$annotation_idx]}
phenotypes=("A3_ratio" "A5_ratio" "AL_ratio" "AF_ratio" "MX_ratio" "RI_ratio" "SE_ratio")
for phenotype in ${phenotypes[@]}
do
python ${py_script} --cohort ${cohort} --annotation ${annotation} --phenotype ${phenotype}
done

covert tensorqtl output to query format
chrs=( {1..22} )
phenotypes=("A3_ratio" "A5_ratio" "AL_ratio" "AF_ratio" "MX_ratio" "RI_ratio" "SE_ratio")
task_id=${SLURM_ARRAY_TASK_ID}
chr_idx=$(( task_id / 7 ))
pheno_idx=$(( task_id % 7 ))
chr=${chrs[$chr_idx]}
phenotype=${phenotypes[$pheno_idx]}
echo "Task ID: ${task_id}, chr: ${chr}, phenotype: ${phenotype}"
py_script=03_metabrain_cohort_data/script/conver_tensorqtl_res_to_SMRquery_refine.py
python ${py_script} \
    --chr ${chr} \
    --phenotype ${phenotype}

# split bed file 
cohorts=("HBCC" "LIBD" "BM10" "BM22" "BM36" "BM44" "ROSMAP" "BrainGVEx" "CMC" "Mayo")
cohort=${cohorts[$SLURM_ARRAY_TASK_ID-1]}
Residual_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/${cohort}/splicing_phenotype/phenotype_Residual
split_script=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/script/split_bedfile.py
python ${split_script} --Residual ${Residual_dir}