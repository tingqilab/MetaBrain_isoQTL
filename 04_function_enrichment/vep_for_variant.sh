#!/bin/bash
#SBATCH -J vep
#SBATCH -N 1
#SBATCH -n 1  
#SBATCH --cpus-per-task=4
#SBATCH -p qitinglab
#SBATCH -e vep_%A_%a.err
#SBATCH -o vep_%A_%a.out
#SBATCH -A qitinglab-grp
#SBATCH --mem=10G
#SBATCH --array=1-22
#SBATCH --nodelist=qiting-svr2

source ~/.bashrc
conda activate vep_env

## vep annotation for denovo mutation
# denovo_autosome='/qitinglab/sharedata/De_novo_mutation_DB/denovo_clean/denovo_autosome.nochr.vcf.gz'
# gtf='/qitinglab/sharedata/genome_annotation_reference/gencode_for_vep.gtf.gz'
chr=$SLURM_ARRAY_TASK_ID
output_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/function_annotation/vep_annot
DNM_vcf=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/all_variant/vcf_file/chr${chr}.vcf
Metabrain_Augment_gtf='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/01_construct_transcript_panel/transcript_panel/SQANTI_QC_full_length/final_panel/post_qc/Metabrain_Augment_final_transcript_panel_CDS.sort.gtf.gz'
ref_fa='/qitinglab/sharedata/genome_annotation_reference/GRCh38.primary_assembly.genome.fa'
chr_synonyms='/qitinglab/sharedata/homo_sapiens/115_GRCh38/chr_synonyms.txt'
cache='/qitinglab/sharedata/homo_sapiens/115_GRCh38'
gencode_gtf=/qitinglab/sharedata/genome_annotation_reference/gencode_for_vep.gtf.gz
vep -i ${DNM_vcf} --dir_cache ${cache} --gtf ${Metabrain_Augment_gtf} --fasta ${ref_fa} -synonyms ${chr_synonyms} --fork 4 \
-o ${output_dir}/chr${chr}.txt --everything