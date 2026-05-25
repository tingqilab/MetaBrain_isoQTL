#!/bin/sh
#SBATCH -J gwas_qtl_coloc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH -p qitinglab
#SBATCH -e gwas_qtl_coloc_%A_%a.err
#SBATCH -o gwas_qtl_coloc_%A_%a.out
#SBATCH -A qitinglab-grp
#SBATCH --mem=10G
#SBATCH --array=0-199

source ~/.bashrc
conda activate jupyter
R_code=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_unique_isoQTL_genetic_control/new_script/COLOC_qlt_gwas.R
out_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/result/GWAS_QTL_coloc_res
mkdir -p ${out_dir}
mkdir -p ${out_dir}/result
mkdir -p ${out_dir}/gene
input_list_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/gwas_qtl_coloc_input_list
CHUNK_ID=$(printf "%03d" $SLURM_ARRAY_TASK_ID)
input_list=${input_list_dir}/gwas_qtl_coloc_list_part_${CHUNK_ID}
for input in $(cat $input_list); do
    echo "[$(date)] Running coloc for input ${input}"
    basename=$(basename $input coloc_input.txt)
    Rscript $R_code --qtl_gwas ${input} --out_results ${out_dir}/result/${basename}coloc_res.txt --out_genes ${out_dir}/gene/${basename}gene.txt --pp4 0.8
done