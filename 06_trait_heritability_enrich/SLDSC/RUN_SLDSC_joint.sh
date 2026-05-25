#!/bin/bash
#SBATCH -J SLDSC_joint
#SBATCH -p qitinglab
#SBATCH -o SLDSC_joint_%A_%a.out
#SBATCH -e SLDSC__joint_%A_%a.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A qitinglab-grp
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --array=1-16   

#### 环境与资源
source ~/.bashrc
conda activate LDSC

baseline_v1=/qitinglab/sharedata/LDSC_resource/GRCh38/baseline_v1.2/baseline.
baseline_v2=/qitinglab/sharedata/LDSC_resource/GRCh38/baselineLD_v2.2/baselineLD.
weights=/qitinglab/sharedata/LDSC_resource/GRCh38/weights/weights.hm3_noMHC.
LD_1000G_EUR=/qitinglab/sharedata/LDSC_resource/GRCh38/plink_files/1000G.EUR.hg38.
hm3_no_MHC_snp=/qitinglab/sharedata/LDSC_resource/GRCh38/hm3_no_MHC.list.txt
LDSC_DIRT=/qitinglab/software/ldsc_software

GWAS_DIRT=/qitinglab/sharedata/GWAS_summary/brain_disease_gwas_newest_EUR/sumstats
joint_res_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_trait_hertability_enrichment/SLDSC/SLDSC_res
annot_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_trait_hertability_enrichment/SLDSC/obs_annot


PHENO=("" "AD_2022_NatGenet" "ADHD2022_NatGenet" "age_at_menarche" \
  "ALS2021" "Brain_volume_adjHeight_2020_NC" "EA5_NatHumBehav_2024_MAMA" \
  "Insomnia" "IQ_NG_2018" "neuroticism_new" "PD_Nalls_FinnGen_META" \
  "PGC3_SCZ_wave3.european.autosome.public.v3" "PGC_ASD_2019_NatGenet" \
  "PGC_BIP2024_eur_no23andMe" "PGC_MDD_2025_no23andMe_eur" \
  "PGC_PTSD_2024NG_eur" "SmokingInitiation" )
TRAIT_NAME=("" "AD" "ADHD" "AAM" "ALS" "BV" "EA" "Insomnia" "IQ" \
  "Neuroticism" "PD" "SCZ" "ASD" "BIP" "MDD" "PTSD" "SmokingInit")


### joint QTL SLDSC ###
trait=${PHENO[${SLURM_ARRAY_TASK_ID}]}
trait_name=${TRAIT_NAME[${SLURM_ARRAY_TASK_ID}]}
trait_gwas=${GWAS_DIRT}/${trait}.sumstats.gz
python ${LDSC_DIRT}/ldsc.py --h2 ${trait_gwas} \
        --ref-ld-chr ${baseline_v1},${annot_dir}/known_isoQTL_obs.,${annot_dir}/novel_isoQTL_obs.,${annot_dir}/eQTL_obs. \
        --w-ld-chr ${weights} --overlap-annot --print-cov --print-coefficients --print-delete-vals --frqfile-chr ${LD_1000G_EUR} \
        --out ${joint_res_dir}/${trait_name}_known_novel_isoQTL_eQTL

python ${LDSC_DIRT}/ldsc.py --h2 ${trait_gwas} \
        --ref-ld-chr ${baseline_v1},${annot_dir}/known_irQTL_obs.,${annot_dir}/novel_irQTL_obs.,${annot_dir}/eQTL_obs. \
        --w-ld-chr ${weights} --overlap-annot --print-cov --print-coefficients --print-delete-vals --frqfile-chr ${LD_1000G_EUR} \
        --out ${joint_res_dir}/${trait_name}_known_novel_irQTL_eQTL

# python ${LDSC_DIRT}/ldsc.py --h2 ${trait_gwas} \
#         --ref-ld-chr ${baseline_v1},${annot_dir}/known_irQTL_obs.,${annot_dir}/novel_irQTL_obs. \
#         --w-ld-chr ${weights} --overlap-annot --print-cov --print-coefficients --print-delete-vals --frqfile-chr ${LD_1000G_EUR} \
#         --out ${joint_res_dir}/${trait_name}_known_novel_irQTL

### irQTL,isoQTL,eQTL,sQTL
# python ${LDSC_DIRT}/ldsc.py --h2 ${trait_gwas} \
#         --ref-ld-chr ${baseline_v1},${annot_dir}/eQTL_obs. \
#         --w-ld-chr ${weights} --overlap-annot --print-cov --print-coefficients --print-delete-vals --frqfile-chr ${LD_1000G_EUR} \
#         --out ${joint_res_dir}/${trait_name}_eQTL
# ### irQTL,isoQTL,eQTL

# python ${LDSC_DIRT}/ldsc.py --h2 ${trait_gwas} \
#         --ref-ld-chr ${baseline_v1},${annot_dir}/sQTL_obs. \
#         --w-ld-chr ${weights} --overlap-annot --print-cov --print-coefficients --print-delete-vals --frqfile-chr ${LD_1000G_EUR} \
#         --out ${joint_res_dir}/${trait_name}_sQTL

# python ${LDSC_DIRT}/ldsc.py --h2 ${trait_gwas} \
#         --ref-ld-chr ${baseline_v1},${annot_dir}/irQTL_obs. \
#         --w-ld-chr ${weights} --overlap-annot --print-cov --print-coefficients --print-delete-vals --frqfile-chr ${LD_1000G_EUR} \
#         --out ${joint_res_dir}/${trait_name}_irQTL

### isoQTL,eQTL,sQTL


# python ${LDSC_DIRT}/ldsc.py --h2 ${trait_gwas} \
#         --ref-ld-chr ${baseline_v1},${annot_dir}/irQTL_obs.,${annot_dir}/eQTL_obs.,${annot_dir}/sQTL_obs. \
#         --w-ld-chr ${weights} --overlap-annot --print-cov --print-coefficients --print-delete-vals --frqfile-chr ${LD_1000G_EUR} \
#         --out ${joint_res_dir}/${trait_name}_irQTL_eQTL_sQTL
# ### irQTL,eQTL,sQTL
# python ${LDSC_DIRT}/ldsc.py --h2 ${trait_gwas} \
#         --ref-ld-chr ${baseline_v1},${annot_dir}/known_irQTL_obs.,${annot_dir}/known_isoQTL_obs.,${annot_dir}/novel_irQTL_obs.,${annot_dir}/novel_isoQTL_obs. \
#         --w-ld-chr ${weights} --overlap-annot --print-cov --print-coefficients --print-delete-vals --frqfile-chr ${LD_1000G_EUR} \
#         --out ${joint_res_dir}/${trait_name}_known_novel_irQTL_isoQTL
# ## known_irQTL,known_isoQTL,novel_irQTL,novel_isoQTL


