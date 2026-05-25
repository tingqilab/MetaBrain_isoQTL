#!/bin/bash

Chresidualdir_list=bed_dir_list.txt
Chresidualdir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $Chresidualdir_list)
cohort=$(basename "$Chresidualdir" | cut -d'_' -f1)
echo "Task $SLURM_ARRAY_TASK_ID"
echo "Chresidualdir: $Chresidualdir"
echo "cohort: $cohort"

tensorqtl_script=eQTL_analysis_permutation_test.py
chr=$SLURM_ARRAY_TASK_ID
Chresidualdir='CMC_Gencode_transcript_tier_PC65'
echo "Processing chr$chr..."
python3 $tensorqtl_script \
    --cohort 'CMC' \
    --Chresidualdir $Chresidualdir \
    --chr $chr --outdir "tensorqtl_out_permutation"