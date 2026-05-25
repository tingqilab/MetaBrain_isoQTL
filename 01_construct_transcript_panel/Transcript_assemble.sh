###############################
# Transcript assembly methods #
###############################

source ~/.bashrc
conda activate flair

bam=input.sorted.bam
outdir=output_directory

################################
# Bambu transcript assembly #
################################

Rscript Bambu.R \
    --bam ${bam} \
    --output_dir ${outdir}

################################
# StringTie3 transcript assembly #
################################

annotation=gencode.v47.annotation.gtf

stringtie \
    -L \
    -G ${annotation} \
    -p 12 \
    -o ${outdir}/stringtie3.out.gtf \
    ${bam}

###############################
# IsoQuant transcript assembly #
###############################

reference=GRCh38.primary_assembly.genome.fa
annotation_db=gencode.v47.annotation.db

isoquant.py \
    --reference ${reference} \
    --bam ${bam} \
    --genedb ${annotation_db} \
    --data_type nanopore or pacbio_ccs \
    --threads 10 \
    -o ${outdir}/IsoQuant \
    --prefix sample \
    --process_only_chr \
        chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
        chr20 chr21 chr22 chrX chrY chrM