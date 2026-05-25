############################
# Oxford Nanopore RNA-seq #
############################
threads=24
reference=GRCh38.primary_assembly.genome.fa
# FASTQ statistics
seqkit stats --all input.fastq > fastq_stats.txt -j ${threads}
# Alignment
minimap2 \
    -ax splice \
    -uf \
    -k14 \
    --secondary=no \
    -t ${threads} \
    ${reference} \
    input.fastq | \
    samtools sort -@ ${threads} -o aligned.sorted.bam

# BAM indexing and statistics
samtools index aligned.sorted.bam
samtools stats \
    -@ ${threads} \
    aligned.sorted.bam > bam_stats.txt

########################
# PacBio Iso-Seq data #
########################
threads=48
reference=GRCh38.primary_assembly.genome.fa

# FASTQ statistics
seqkit stats --all input.fastq > fastq_stats.txt
# Alignment
minimap2 \
    -ax splice:hq \
    -C5 \
    --secondary=no \
    -t ${threads} \
    ${reference} \
    input.fastq | \
    samtools sort -@ ${threads} -o aligned.sorted.bam

# BAM indexing and statistics
samtools index aligned.sorted.bam
samtools stats \
    -@ ${threads} \
    aligned.sorted.bam > bam_stats.txt