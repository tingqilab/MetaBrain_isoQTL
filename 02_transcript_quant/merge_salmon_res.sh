## Merge Salmon transcript quantification results
# This script aggregates Salmon quantification outputs from multiple
# transcriptome panels using tximport.
#
# Expected directory structure:
#
# ${work_dir}/${cohort_name}/
# ├── Gencode/
# │   ├── sample1/
# │   ├── sample2/
# │   └── ...
# ├── Metabrain/
# │   ├── sample1/
# │   ├── sample2/
# │   └── ...
# └── Metabrain_Augment/
#     ├── sample1/
#     ├── sample2/
#     └── ...
#
# The script:
#   1. Finds all Salmon quant.sf files
#   2. Aggregates transcript quantifications using tximport
#   3. Generates merged abundance matrices
#   4. Compresses intermediate result files

#!/bin/sh

### input paths
work_dir=/path/to/work_dir
tx2gene_path=/path/to/tx2gene.txt
aggregate_Rscript=/path/to/taximport.R

### cohort
cohorts=('ROSMAP')
cohort_name=${cohorts[0]}

cd ${work_dir}/${cohort_name}
cohort_dir=${work_dir}/${cohort_name}

### collect Salmon quantification files
find "${cohort_dir}/Gencode" \
    -maxdepth 2 \
    -name "quant.sf" \
    > ./Gencode_${cohort_name}_mergequant_file.txt

find "${cohort_dir}/Metabrain" \
    -maxdepth 2 \
    -name "quant.sf" \
    > ./Metabrain_${cohort_name}_mergequant_file.txt

find "${cohort_dir}/Metabrain_Augment" \
    -maxdepth 2 \
    -name "quant.sf" \
    > ./Metabrain_Augment_${cohort_name}_mergequant_file.txt

### output directory
out_put_dir=${cohort_dir}

source ~/.bashrc
conda activate jupyter

### merge transcript quantifications
Rscript ${aggregate_Rscript} \
    ${tx2gene_path} \
    ./Gencode_${cohort_name}_mergequant_file.txt \
    ${out_put_dir} \
    ${cohort_name}_Gencode

Rscript ${aggregate_Rscript} \
    ${tx2gene_path} \
    ./Metabrain_${cohort_name}_mergequant_file.txt \
    ${out_put_dir} \
    ${cohort_name}_Metabrain

Rscript ${aggregate_Rscript} \
    ${tx2gene_path} \
    ./Metabrain_Augment_${cohort_name}_mergequant_file.txt \
    ${out_put_dir} \
    ${cohort_name}_Metabrain_Augment

### clean temporary file lists
rm ./Gencode_${cohort_name}_mergequant_file.txt
rm ./Metabrain_${cohort_name}_mergequant_file.txt
rm ./Metabrain_Augment_${cohort_name}_mergequant_file.txt

### archive merged result files
mkdir -p ./merge_result

mv *.txt ./merge_result

tar -czvf merge_result.tar.gz merge_result

rm -r merge_result


