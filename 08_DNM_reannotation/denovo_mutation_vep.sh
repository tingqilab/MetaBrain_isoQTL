#!/bin/sh
#SBATCH -J vep_annotation
#SBATCH -p qitinglab
#SBATCH --output=vep_annotation.out
#SBATCH --error=vep_annotation.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A qitinglab-grp
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --nodelist=qiting-svr1

source ~/.bashrc
conda activate vep_env

## vep annotation for denovo mutation
# denovo_autosome='/qitinglab/sharedata/De_novo_mutation_DB/denovo_clean/denovo_autosome.nochr.vcf.gz'
# gtf='/qitinglab/sharedata/genome_annotation_reference/gencode_for_vep.gtf.gz'

DNM_vcf=/qitinglab/sharedata/De_novo_mutation_DB/new_try/Gene4Denovo2/DNM.vcf.gz
clinvar=/qitinglab/sharedata/De_novo_mutation_DB/new_try/clinvar/clinvar.vcf.gz
Metabrain_Augment_gtf='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/01_construct_transcript_panel/transcript_panel/SQANTI_QC_full_length/final_panel/post_qc/Metabrain_Augment_final_transcript_panel_CDS_update.sort.gtf.gz'
ref_fa='/qitinglab/sharedata/genome_annotation_reference/GRCh38.primary_assembly.genome.fa'
chr_synonyms='/qitinglab/sharedata/homo_sapiens/115_GRCh38/chr_synonyms.txt'
cache='/qitinglab/sharedata/homo_sapiens/115_GRCh38'
gencode_gtf=/qitinglab/sharedata/genome_annotation_reference/gencode_for_vep.gtf.gz
vep -i ${DNM_vcf} --dir_cache ${cache} --gtf ${Metabrain_Augment_gtf} --fasta ${ref_fa} -synonyms ${chr_synonyms} --fork 24 \
-o /qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/01_construct_transcript_panel/transcript_final_panel_analysis/re_annot_DNM/DNM_Metabrain_Augment_Vepout_update.txt \
--everything


##  --fork 8 : vep multi threads parameter
## vep annotation for control population variant gnomAD DB
# gnomAD_dir='/qitinglab/sharedata/genomeAD/genomeAD_assess_selection'
# vep_result_dir='/qitinglab/sharedata/genomeAD/genomeAD_assess_selection/vep_result'