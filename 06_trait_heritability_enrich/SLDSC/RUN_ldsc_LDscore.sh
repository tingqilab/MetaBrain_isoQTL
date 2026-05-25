#!/bin/bash
#SBATCH -J compute_LD_annot
#SBATCH -p qitinglab
#SBATCH -o computeannotLDscore_%A_%a.out
#SBATCH -e computeannotLDscore_%A_%a.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A qitinglab-grp
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --nodelist=qiting-svr2
#SBATCH --array=0-200%100

### 0-17599
#### LDSC_resoure and software environment
baseline_v1=/qitinglab/sharedata/LDSC_resource/GRCh38/baseline_v1.2/baseline.
baseline_v2=/qitinglab/sharedata/LDSC_resource/GRCh38/baselineLD_v2.2/baselineLD.
weights=/qitinglab/sharedata/LDSC_resource/GRCh38/weights/weights.hm3_noMHC.
LD_1000G_EUR=/qitinglab/sharedata/LDSC_resource/GRCh38/plink_files/1000G.EUR.hg38
hm3_no_MHC_snp=/qitinglab/sharedata/LDSC_resource/GRCh38/hm3_no_MHC.list.txt
LDSC_DIRT=/qitinglab/software/ldsc_software
source ~/.bashrc
conda activate LDSC
#### LDSC_resoure and software environment

### caculate the LD for obs snp of topQTL
# obs_annot_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_trait_hertability_enrichment/SLDSC/obs_annot
# types=("novel_isoQTL" "known_isoQTL" "isoQTL" "novel_irQTL" "known_irQTL" "irQTL" "eQTL" "sQTL")
# chrs=($(seq 1 22))
# type_index=$(( (SLURM_ARRAY_TASK_ID - 1) / 22 ))
# chr_index=$(( (SLURM_ARRAY_TASK_ID - 1) % 22 ))
# type=${types[$type_index]}
# chr=${chrs[$chr_index]}
# echo "[$(date)] Running caculate LD score type=${type}, chr=${chr}"
# if [ -f "${obs_annot_dir}/${type}_obs.${chr}.l2.ldscore.gz" ]; then
#     echo "⚠️  Output file already exists: ${out_prefix}.l2.ldscore.gz — skipping."
#     exit 0
# fi

# python ${LDSC_DIRT}/ldsc.py --l2 \
#     --bfile ${LD_1000G_EUR}.${chr} \
#     --annot ${obs_annot_dir}/${type}_obs.${chr}.annot.gz \
#     --out ${obs_annot_dir}/${type}_obs.${chr} \
#     --ld-wind-cm 1 \
#     --print-snps ${hm3_no_MHC_snp}

### caculate the LD for null snp of topQTL
null_annot_dir=/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_trait_hertability_enrichment/SLDSC/null_annot
# ====== 参数定义 ======
types=("eQTL" "sQTL" "isoQTL" "irQTL" "known_isoQTL" "novel_isoQTL" "known_irQTL" "novel_irQTL")
chrs=($(seq 1 22))
reps=($(seq 1 100))

total_types=${#types[@]}    # 8
total_chrs=${#chrs[@]}      # 22
total_reps=${#reps[@]}      # 100

# ====== 当前任务编号 ======
task_id=$SLURM_ARRAY_TASK_ID

# ====== 控制顺序为 rep -> type -> chr ======
rep_index=$(( task_id / (total_types * total_chrs) ))
type_index=$(( (task_id / total_chrs) % total_types ))
chr_index=$(( task_id % total_chrs ))

rep=${reps[$rep_index]}
type=${types[$type_index]}
chr=${chrs[$chr_index]}

echo "[$(date '+%F %T')] Task $SLURM_ARRAY_TASK_ID -> rep=$rep, type=$type, chr=$chr"

# ====== 输出路径 ======
out_prefix=${null_annot_dir}/${type}_rep${rep}.${chr}
out_file=${out_prefix}.l2.ldscore.gz

# ====== 跳过已完成文件 ======
if [ -f "$out_file" ]; then
    echo "Output exists, skipping: $out_file"
    exit 0
fi

# ====== 运行LDSC ======
python ${LDSC_DIRT}/ldsc.py --l2 \
    --bfile ${LD_1000G_EUR}.${chr} \
    --annot ${null_annot_dir}/${type}_rep${rep}.${chr}.annot.gz \
    --out ${out_prefix} \
    --ld-wind-cm 1 \
    --print-snps ${hm3_no_MHC_snp}

echo "[$(date '+%F %T')] Finished -> ${out_prefix}"
