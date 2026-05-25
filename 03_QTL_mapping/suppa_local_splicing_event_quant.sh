#!/bin/sh
#SBATCH -J suppa
#SBATCH -p qitinglab
#SBATCH -e suppa_%A_%a.err
#SBATCH -o suppa_%A_%a.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A qitinglab-grp
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --array=1-10

source ~/.bashrc
conda activate flair
cohorts=("HBCC" "LIBD" "BM10" "BM22" "BM36" "BM44" "ROSMAP" "BrainGVEx" "CMC" "Mayo")
cohort=${cohorts[$SLURM_ARRAY_TASK_ID-1]}
suppa_py=/qitinglab/software/SUPPA-2.4/suppa.py

input_dir=03_metabrain_cohort_data
transcript_expression_dir=${input_dir}/${cohort}/phenotype/merge_result
Gencode_expression=${transcript_expression_dir}/${cohort}_Gencode_transcript_abundance.txt
Metabrain_expression=${transcript_expression_dir}/${cohort}_Metabrain_transcript_abundance.txt
Metabrain_Augment_expression=${transcript_expression_dir}/${cohort}_Metabrain_Augment_transcript_abundance.txt
out_dir=${input_dir}/${cohort}/splicing_phenotype
mkdir -p ${out_dir}
mkdir -p ${out_dir}/tmp
cd ${out_dir}
awk 'NR==1 {sub("^[^\t]+\t","")}1' ${Gencode_expression} > ./tmp/Gencode_transcript.txt
awk 'NR==1 {sub("^[^\t]+\t","")}1' ${Metabrain_expression} > ./tmp/Metabrain_transcript.txt
awk 'NR==1 {sub("^[^\t]+\t","")}1' ${Metabrain_Augment_expression} > ./tmp/Metabrain_Augment_transcript.txt

Gencode_splice_event_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/01_construct_transcript_panel/transcript_final_panel_analysis/suppa_splice_events/Gencode
Metabrain_splice_event_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/01_construct_transcript_panel/transcript_final_panel_analysis/suppa_splice_events/Metabrain
Metabrain_Augment_splice_event_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/01_construct_transcript_panel/transcript_final_panel_analysis/suppa_splice_events/Metabrain_Augment
rm -r ${out_dir}/suppa_splicing 
mkdir -p ${out_dir}/suppa_splicing
final_out_dir=${out_dir}/suppa_splicing
event_type=('A3' 'A5' 'AF' 'AL' 'MX' 'RI' 'SE')
for event in "${event_type[@]}"
do
echo "$event"
Gencode_event_ioe=${Gencode_splice_event_dir}/Gencode_${event}_strict.ioe
Metabrain_event_ioe=${Metabrain_splice_event_dir}/Metabrain_${event}_strict.ioe
Metabrain_Augment_event_ioe=${Metabrain_Augment_splice_event_dir}/Metabrain_Augment_${event}_strict.ioe
python ${suppa_py} psiPerEvent -i ${Gencode_event_ioe} -e ./tmp/Gencode_transcript.txt -o ${final_out_dir}/${cohort}_Gencode_${event}
python ${suppa_py} psiPerEvent -i ${Metabrain_event_ioe} -e ./tmp/Metabrain_transcript.txt -o ${final_out_dir}/${cohort}_Metabrain_${event}
python ${suppa_py} psiPerEvent -i ${Metabrain_Augment_event_ioe} -e ./tmp/Metabrain_Augment_transcript.txt -o ${final_out_dir}/${cohort}_Metabrain_Augment_${event}
done
#rm -r ${out_dir}/tmp