#!/bin/sh
#SBATCH -J calculd
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH -p qitinglab
#SBATCH -e calculd_%A_%a.err
#SBATCH -o calculd_%A_%a.out
#SBATCH -A qitinglab-grp
#SBATCH --mem=6G
#SBATCH --array=1-22
#SBATCH --nodelist=qiting-svr1

chr=${SLURM_ARRAY_TASK_ID}
plink=/qitinglab/software/plink
lead_QTL_pair_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_unique_isoQTL_genetic_control/leadQTL_snp_pair
bfile_hg38_1000G=/qitinglab/sharedata/LD_reference/1000G/1000G_Phase3_EUR_plinkfiles/1000G.EUR.QC.
lead_QTL_pair=${lead_QTL_pair_dir}/chr${SLURM_ARRAY_TASK_ID}_isoQTL_isoQTL_ld.txt
out_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_unique_isoQTL_genetic_control/lead_SNP_ld_res
output=${out_dir}/chr${SLURM_ARRAY_TASK_ID}_isoQTL_isoQTL_ld_res.txt
cd /qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_unique_isoQTL_genetic_control
cut -f3,4 ${lead_QTL_pair} | tr '\t' '\n' | sort -u > chr${chr}_snplist.txt

${plink} \
--bfile ${bfile_hg38_1000G}${chr} \
--extract chr${chr}_snplist.txt \
--r2 \
--ld-window-kb 10000 \
--ld-window 999999 \
--ld-window-r2 0 \
--out ${out_dir}/chr${chr}_LD

rm chr${chr}_snplist.txt