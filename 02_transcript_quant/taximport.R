library(tximport)
library(edgeR)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop(
        "Usage:\n",
        "Rscript tximport_gene_transcript_aggregate.R ",
        "<tx2gene.txt> <salmon_quant_list.txt> <out_dir> <out_prefix>\n\n",
        "Example:\n",
        "Rscript tximport_gene_transcript_aggregate.R ",
        "tx2gene.txt Gencode_CMC_quant_file.txt merge_data CMC_Gencode\n"
    )
}

tx2gene_file     <- args[1]
salmon_list_file <- args[2]
out_dir          <- args[3]
out_prefix       <- args[4]


tx2gene <- read.table(
    tx2gene_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

salmon_paths <- readLines(salmon_list_file)
names(salmon_paths) <- basename(dirname(salmon_paths))

# -----------------------------
# 3. tximport
# -----------------------------
txi_gene_level <- tximport(
    salmon_paths,
    type = "salmon",
    tx2gene = tx2gene,
    ignoreAfterBar = TRUE
)

txi_transcript_level <- tximport(
    salmon_paths,
    type = "salmon",
    txOut = TRUE
)

# -----------------------------
#  gene-level
# -----------------------------
gene_counts    <- txi_gene_level$counts
gene_abundance <- txi_gene_level$abundance
gene_length    <- txi_gene_level$length

gene_counts_df <- cbind(
    gene_id = rownames(gene_counts),
    as.data.frame(gene_counts)
)

gene_abundance_df <- cbind(
    gene_id = rownames(gene_abundance),
    as.data.frame(gene_abundance)
)

gene_length_df <- cbind(
    gene_id = rownames(gene_length),
    as.data.frame(gene_length)
)

# -----------------------------
#  transcript-level
# -----------------------------
transcript_counts    <- txi_transcript_level$counts
transcript_abundance <- txi_transcript_level$abundance
transcript_length    <- txi_transcript_level$length

transcript_counts_df <- cbind(
    transcript_id = rownames(transcript_counts),
    as.data.frame(transcript_counts)
)

transcript_abundance_df <- cbind(
    transcript_id = rownames(transcript_abundance),
    as.data.frame(transcript_abundance)
)

transcript_length_df <- cbind(
    transcript_id = rownames(transcript_length),
    as.data.frame(transcript_length)
)


####### correct RTA by resampling transcript counts #############
salmon_paths_parent <- dirname(salmon_paths)
y <- catchSalmon(salmon_paths_parent, verbose=TRUE)
dge <- DGEList(counts=y$counts/y$annotation$Overdispersion, genes=y$annotation)
Scaled_Count<-dge$counts
new_colnames <- basename(salmon_paths_parent)
colnames(Scaled_Count) <- new_colnames
Scaled_Count_df <- as.data.frame(Scaled_Count)
Scaled_Count_df <- cbind(transcript_id = rownames(Scaled_Count_df), Scaled_Count_df)

# -----------------------------
# 6. output
# -----------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.table(
    gene_counts_df,
    file = file.path(out_dir, paste0(out_prefix, "_gene_counts.txt")),
    sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
    gene_abundance_df,
    file = file.path(out_dir, paste0(out_prefix, "_gene_abundance.txt")),
    sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
    gene_length_df,
    file = file.path(out_dir, paste0(out_prefix, "_gene_length.txt")),
    sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
    transcript_counts_df,
    file = file.path(out_dir, paste0(out_prefix, "_transcript_counts.txt")),
    sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
    transcript_abundance_df,
    file = file.path(out_dir, paste0(out_prefix, "_transcript_abundance.txt")),
    sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
    transcript_length_df,
    file = file.path(out_dir, paste0(out_prefix, "_transcript_length.txt")),
    sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(Scaled_Count_df, 
    file = file.path(out_dir, paste0(out_prefix, "_transcript_RTAcorrect_count.txt")), 
    row.names = FALSE,sep = "\t",quote = FALSE
)
