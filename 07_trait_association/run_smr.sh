#!/bin/bash
#SBATCH -J SMR
#SBATCH -p qitinglab
#SBATCH -A qitinglab-grp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --array=1-352%20
#SBATCH --nodelist=qiting-svr1
#SBATCH -e SMR_%A_%a.err
#SBATCH -o SMR_%A_%a.out

smr=/qitinglab/software/smr_Linux
GWAS_DIR="/qitinglab/sharedata/GWAS_summary/brain_disease_gwas_newest_EUR"
REFERENCE="/qitinglab/sharedata/LD_reference/1000G/1000G_Phase3_EUR_plinkfiles/1000G.EUR.QC."
SMR_OUTPUT="/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/SMR/smr_res"
#SMR_OUTPUT_multi='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/SMR/smr_multi_res'
mkdir -p ${SMR_OUTPUT}
#mkdir -p ${SMR_OUTPUT_multi}
QTL_besd_dir='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/Meta_analysis/Meta_QTL/besd'
# format : Metabrain_Augment_A5_ratio_chr20.besd
# annotation list
annotations=("Gencode" "Metabrain_Augment")
# QTL type list
qtl_types=('eQTL' 'sQTL' 'isoQTL' 'irQTL' 'A3_ratio' 'A5_ratio' 'AL_ratio' 'AF_ratio' 'MX_ratio' 'RI_ratio' 'SE_ratio')
# trait GWAS file base name
PHENO=("AD_2022_NatGenet" "ADHD2022_NatGenet" "age_at_menarche" \
  "ALS2021" "Brain_volume_adjHeight_2020_NC" "EA5_NatHumBehav_2024_MAMA" \
  "Insomnia" "IQ_NG_2018" "neuroticism_new" "PD_Nalls_FinnGen_META" \
  "PGC3_SCZ_wave3.european.autosome.public.v3" "PGC_ASD_2019_NatGenet" \
  "PGC_BIP2024_eur_no23andMe" "PGC_MDD_2025_no23andMe_eur" \
  "PGC_PTSD_2024NG_eur" "SmokingInitiation")
trait_name=('AD' 'ADHD' 'AAM' 'ALS' 'BV' 'EA5' 'Ins' 'IQ' 'Neu' 'PD' 'SCZ' 'ASD' 'BIP' 'MDD' 'PTSD' 'Smoking')

############################
# decode index to get annotation, qtl type and trait
############################
ann_n=${#annotations[@]}
qtl_n=${#qtl_types[@]}
pheno_n=${#PHENO[@]}

task_id=$((SLURM_ARRAY_TASK_ID - 1))

ann_index=$(( task_id / (qtl_n * pheno_n) ))
qtl_index=$(( (task_id / pheno_n) % qtl_n ))
pheno_index=$(( task_id % pheno_n ))

annotation=${annotations[$ann_index]}
qtl=${qtl_types[$qtl_index]}
pheno=${PHENO[$pheno_index]}
trait=${trait_name[$pheno_index]}

echo "annotation: $annotation"
echo "qtl: $qtl"
echo "pheno: $pheno"
echo "trait: $trait"

for chr in {1..22}; do
    QTL_FILE="${QTL_besd_dir}/${annotation}_${qtl}_chr${chr}" 
    echo "Running SMR for chromosome $chr"
    SMR_OUT_PREFIX="${SMR_OUTPUT}/${annotation}_${qtl}_${trait}_chr${chr}"
    source ~/.bashrc
    conda activate smr_env
    # run SMR with smr multi
    $smr --bfile ${REFERENCE}${chr} \
         --beqtl-summary $QTL_FILE \
         --gwas-summary ${GWAS_DIR}/${pheno}.txt \
         --maf 0.01 \
         --out ${SMR_OUT_PREFIX} \
         --thread-num 4 --smr-multi    
done
