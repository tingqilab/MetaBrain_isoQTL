#!/bin/sh
#SBATCH -J suppa
#SBATCH -p qitinglab
#SBATCH -e suppa.erro
#SBATCH -o suppa.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A qitinglab-grp
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --nodelist=qiting-svr1
#SBATCH --array=0-2

source ~/.bashrc
conda activate flair

suppa_py=/qitinglab/software/SUPPA-2.4/suppa.py
gtf_panel_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/01_construct_transcript_panel/transcript_final_panel_new
suppa_out_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/01_construct_transcript_panel/transcript_final_panel_analysis/suppa_splice_events
cd /qitinglab/shenlongxin/metabrain/NOVEL_DISCOVERY/final_panel/qc_result_final/post_qc/suppa
sample_list=('Gencode' 'Metabrain' 'Metabrain_Augment')
sample=${sample_list[$SLURM_ARRAY_TASK_ID]}
gtf=${gtf_panel_dir}/${sample}/${sample}.gtf
mkdir -p ${suppa_out_dir}/${sample}
cd ${suppa_out_dir}/${sample} 
python ${suppa_py} generateEvents -i ${gtf} -o ${sample} -f ioe -e RI SE SS MX FL