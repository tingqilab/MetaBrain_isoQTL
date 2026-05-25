#!/bin/sh
#SBATCH -J clumping
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH -p qitinglab
#SBATCH -e clumping_%A_%a.err
#SBATCH -o clumping_%A_%a.out
#SBATCH -A qitinglab-grp
#SBATCH --mem=10G
#SBATCH --nodelist=qiting-svr2
#SBATCH --array=1-16

plink=/qitinglab/software/plink
GWAS_dir=/qitinglab/sharedata/GWAS_summary/brain_disease_gwas_newest_EUR
hg37_out_dir=/qitinglab/sharedata/GWAS_summary/brain_disease_gwas_newest_EUR/clumping_res_hg37_LD
hg38_out_dir=/qitinglab/sharedata/GWAS_summary/brain_disease_gwas_newest_EUR/clumping_res_hg38_LD

mkdir -p ${hg37_out_dir}
mkdir -p ${hg38_out_dir}

bfile_hg38_1000G=/qitinglab/sharedata/LD_reference/1000G/1000G_Phase3_EUR_plinkfiles/1000G.EUR.QC.
bfile_hg19_1000G=/qitinglab/chenhuilin/data/1000G_genotype/1000G_eur_split_QC/1000G.EUR.QC.


gwas_names=("AD_2022_NatGenet" "ADHD2022_NatGenet" "age_at_menarche" \
  "ALS2021" "Brain_volume_adjHeight_2020_NC" "EA5_NatHumBehav_2024_MAMA" \
  "Insomnia" "IQ_NG_2018" "neuroticism_new" "PD_Nalls_FinnGen_META" \
  "PGC3_SCZ_wave3.european.autosome.public.v3" "PGC_ASD_2019_NatGenet" \
  "PGC_BIP2024_eur_no23andMe" "PGC_MDD_2025_no23andMe_eur" \
  "PGC_PTSD_2024NG_eur" "SmokingInitiation")


# 获取当前任务对应的性状
index=$(($SLURM_ARRAY_TASK_ID - 1))
trait=${gwas_names[$index]}
gwas_sumstats=${GWAS_dir}/${trait}.txt
echo "Processing trait: $trait (Array task: $SLURM_ARRAY_TASK_ID)"
# 为每个性状创建输出目录
mkdir -p ${hg37_out_dir}/${trait}
mkdir -p ${hg38_out_dir}/${trait}
# To ensure stringent identification of independent signals, we applied the
# following parameters: genome-wide significance thresholds for both index and secondary variants
# 对每个染色体运行clumping
for chr in {1..22}
do
${plink} --bfile ${bfile_hg38_1000G}${chr} \
  --clump ${gwas_sumstats} \
  --clump-p1 5e-8 \
  --clump-p2 5e-8 \
  --clump-r2 0.05 \
  --clump-kb 1000 \
  --clump-field P --clump-snp-field SNP \
  --out ${hg38_out_dir}/${trait}/${trait}_chr${chr}

${plink} --bfile ${bfile_hg19_1000G}${chr} \
  --clump ${gwas_sumstats} \
  --clump-p1 5e-8 \
  --clump-p2 5e-8 \
  --clump-r2 0.05 \
  --clump-kb 1000 \
  --clump-field P --clump-snp-field SNP \
  --out ${hg37_out_dir}/${trait}/${trait}_chr${chr}
done

# 合并所有染色体的结果
cd ${hg37_out_dir}/${trait}
awk 'FNR==1 && NR!=1{next} 1' ${trait}_chr*.clumped | grep -v '^[[:space:]]*$' > ${trait}_all_chr.clumped
echo "Completed processing for trait: $trait"

# 合并所有染色体的结果
cd ${hg38_out_dir}/${trait}
awk 'FNR==1 && NR!=1{next} 1' ${trait}_chr*.clumped | grep -v '^[[:space:]]*$' > ${trait}_all_chr.clumped
echo "Completed processing for trait: $trait"