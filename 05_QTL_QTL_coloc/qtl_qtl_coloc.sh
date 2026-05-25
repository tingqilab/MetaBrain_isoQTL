#!/bin/sh
#SBATCH -J qtl_qtl_coloc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH -p qitinglab
#SBATCH -e qtl_qtl_coloc_%A_%a.err
#SBATCH -o qtl_qtl_coloc_%A_%a.out
#SBATCH -A qitinglab-grp
#SBATCH --mem=10G
#SBATCH --array=1-272

source ~/.bashrc
conda activate jupyter
# R_code=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_unique_isoQTL_genetic_control/new_script/COLOC_isoqtl_isoqtl.R
# out_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/result/isoQTL_isoQTL_coloc_res
# input_list_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/qtl_qtl_coloc_input_list
# CHUNK_ID=$(printf "%03d" $SLURM_ARRAY_TASK_ID)
# input_list=${input_list_dir}/qtl_qtl_coloc_part_${CHUNK_ID}
# for input in $(cat $input_list); do
#     echo "[$(date)] Running coloc for input ${input}"
#     basename=$(basename $input _sig_isoQTL_summary.txt)
#     Rscript $R_code --qtl_qtl ${input} --out_results ${out_dir}/result/${basename}_coloc_res.txt --out_genes ${out_dir}/gene/${basename}.txt --pp4 0.8
# done
## make isoQTL & xQTL COLOC input file
# python /qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_unique_isoQTL_genetic_control/new_script/process_isoQTL_xQTL_input.py ${SLURM_ARRAY_TASK_ID}

## make isoQTL & epiQTL COLOC input file
# python /qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_unique_isoQTL_genetic_control/new_script/process_isoQTL_epiQTL_input.py ${SLURM_ARRAY_TASK_ID}

# isoQTL_xQTL_coloc_R=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_unique_isoQTL_genetic_control/new_script/COLOC_isoqtl_xqtl.R
# isoQTL_xQTL_coloc_res_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/isoQTL_xQTL_input/coloc_res1
# isoQTL_xQTL_sig_coloc_res_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/isoQTL_xQTL_input/sig_coloc_res1
# input_list_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/isoQTL_xQTL_input/split_file_part
# CHUNK_ID=$(printf "%03d" $SLURM_ARRAY_TASK_ID)
# input_list=${input_list_dir}/part${CHUNK_ID}
# for input in $(cat $input_list); do
#     echo "[$(date)] Running coloc for input ${input}"
#     basename=$(basename $input _summary.txt)
#     Rscript $isoQTL_xQTL_coloc_R --qtl_qtl ${input} --out_results ${isoQTL_xQTL_coloc_res_dir}/${basename}_coloc_res.txt --out_genes ${isoQTL_xQTL_sig_coloc_res_dir}/${basename}_sig_coloc_res.txt --pp4 0.8
# done


# chr=$SLURM_ARRAY_TASK_ID
# python /qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_unique_isoQTL_genetic_control/new_script/split_isoQTL_epiQTL_coloc_input.py ${chr}


# -------------------------------
# 获取对应行的后缀
# -------------------------------
index_file=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/isoQTL_epiQTL_input/pair_Probe_index.txt
index=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${index_file})
echo "Processing $index"
work_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/isoQTL_epiQTL_input
Probe_dir=${work_dir}/split_pair_Probe
summary_dir=${work_dir}/split_overlap_summary_data
res_dir=${work_dir}/coloc_res
isoQTL_epiQTL_coloc_R=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_unique_isoQTL_genetic_control/new_script/COLOC_isoqtl_epiqtl.R
Probe_pair_file=${Probe_dir}/Probe_pair_${index}
isoQTL=${summary_dir}/isoQTL_summary_${index}
epiQTL=${summary_dir}/epiQTL_summary_${index}
output=${res_dir}/coloc_res_${index}
Rscript ${isoQTL_epiQTL_coloc_R} -o ${output} -p ${Probe_pair_file} -i ${isoQTL} -e ${epiQTL}
