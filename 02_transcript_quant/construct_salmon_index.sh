## Salmon transcriptome index construction workflow
# This script:
#   1. Generates a decoy list from the reference genome
#   2. Converts transcript annotations (GTF) to transcript FASTA
#   3. Combines transcriptome and genome FASTA files
#   4. Builds a Salmon decoy-aware transcriptome index

#!/bin/sh

### construct decoys.txt
grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | \
    cut -d " " -f 1 > decoys.txt

sed -i.bak -e 's/>//g' decoys.txt

### paths
work_dir=/path/to/work_dir

sample_list=('Metabrain' 'Metabrain_Augment')
sample=${sample_list[$SLURM_ARRAY_TASK_ID]}

gtf=${work_dir}/${sample}/${sample}.gtf
ref_fa=/path/to/reference_genome.fa
decoys=/path/to/decoys.txt

### convert GTF to transcript FASTA
source ~/.bashrc
conda activate gffread

gffread ${gtf} \
    -g ${ref_fa} \
    -w ${work_dir}/${sample}/${sample}.fa

gzip ${work_dir}/${sample}/${sample}.fa

### build Salmon index
conda deactivate
conda activate flair

transcript_fa_gz=${work_dir}/${sample}/${sample}.fa.gz
genome_fa_gz=/path/to/reference_genome.fa.gz

cat ${transcript_fa_gz} ${genome_fa_gz} > \
    ${work_dir}/${sample}/${sample}_combined.fa.gz

salmon index \
    -t ${work_dir}/${sample}/${sample}_combined.fa.gz \
    -i ${work_dir}/${sample}/${sample}_index \
    --decoys ${decoys} \
    -k 31 \
    -p 24
