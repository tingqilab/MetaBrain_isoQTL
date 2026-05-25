## Salmon transcript quantification workflow
# This script performs paired-end RNA-seq transcript quantification
# using Salmon with multiple transcriptome indices.
#
# Features:
#   1. Quantification against multiple transcriptome panels
#   2. Sequence-specific and GC bias correction
#   3. Selective alignment with validateMappings
#   4. Bootstrap-based abundance uncertainty estimation
#
# Recommended resources:
#   - 24 CPU threads
#   - 50 GB RAM

### transcriptome indices
final_panel=/path/to/final_panel

Gencode_index=${final_panel}/Gencode/Gencode_index
Metabrain_Augment_index=${final_panel}/Metabrain_Augment/Metabrain_Augment_index

index_list=(
    ${Metabrain_Augment_index}
    ${Gencode_index}
)

index=${index_list[$SLURM_ARRAY_TASK_ID]}
basename=$(basename ${index} _index)

### paired-end FASTQ files
fq_R1=/path/to/sample_R1.fq.gz
fq_R2=/path/to/sample_R2.fq.gz

### output directory
cd /path/to/output_dir

### run Salmon quantification
salmon quant \
    -i ${index} \
    -l A \
    -1 ${fq_R1} \
    -2 ${fq_R2} \
    -o ./${basename} \
    --seqBias \
    --gcBias \
    --validateMappings \
    --numBootstraps 24 \
    -p 24

### fragment length prior:
# mean = 250
# standard deviation = 25