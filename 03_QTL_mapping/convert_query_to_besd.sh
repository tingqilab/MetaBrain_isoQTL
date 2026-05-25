#!/bin/bash

smr=smr_Linux
input_dir=cohort_query_format
output_dir=cohort_besd_format
annotations=("Metabrain_Augment" "Gencode") # 3 annotation
cohorts=("ROSMAP" "LIBD" "BM10" "BM22" "BM36" "BM44" "HBCC" "BrainGVEx" "CMC" "Mayo") # 10 cohort 
QTL_types=("A3_ratio" "A5_ratio" "AL_ratio" "AF_ratio" "MX_ratio" "RI_ratio" "SE_ratio" "transcript_count" "gene_count" "transcript_ratio") 
chrs=({1..22}) # chr 1-22 22 chr

task_id=$SLURM_ARRAY_TASK_ID
i_ann=$(( task_id / (10 * 7 * 22) ))
rem1=$(( task_id % (10 * 7 * 22) ))
i_coh=$(( rem1 / (7 * 22) ))
rem2=$(( rem1 % (7 * 22) ))
i_qtl=$(( rem2 / 22 ))
i_chr=$(( rem2 % 22 ))


annotation=${annotations[$i_ann]}
cohort=${cohorts[$i_coh]}
qtl=${QTL_types[$i_qtl]}
chr=${chrs[$i_chr]}
echo "Running task: annotation=$annotation, cohort=$cohort, QTL_type=$qtl, chr=$chr"
qfile=${input_dir}/${annotation}/${cohort}_${qtl}_chr${chr}.txt
besd_file=${output_dir}/${annotation}/${cohort}_${qtl}_chr${chr}
case "$cohort" in
    ROSMAP)   n=832 ;;
    CMC)      n=469 ;;
    BrainGVEx) n=384 ;;
    LIBD)     n=128 ;;
    HBCC)     n=129 ;;
    Mayo)     n=288 ;;
    BM10)     n=183 ;;
    BM22)     n=158 ;;
    BM36)     n=138 ;;
    BM44)     n=156 ;;
    *)
        echo "Error: unknown cohort $cohort"
        exit 1
        ;;
esac
${smr} --qfile ${qfile} --make-besd --out ${besd_file} --add-n ${n}


