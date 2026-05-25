#!/bin/bash

Chresidualdir_list=bed_dir_list.txt
Chresidualdir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $Chresidualdir_list)
cohort=$(basename "$Chresidualdir" | cut -d'_' -f1)

echo "Task $SLURM_ARRAY_TASK_ID"
echo "Chresidualdir: $Chresidualdir"
echo "cohort: $cohort"

# ========================
# 循环处理所有 chr
# ========================
tensorqtl_script=eQTL_analysis_only_nominal.py
for chr in {1..22}; do
    echo "Processing chr$chr..."
    python3 $tensorqtl_script \
        --cohort $cohort \
        --Chresidualdir $Chresidualdir \
        --chr $chr --outdir "tensorqtl_out"
done


