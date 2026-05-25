suppressPackageStartupMessages({
  library(bambu)
  library(optparse)
})
option_list <- list(
  make_option(c("--bam"), type = "character", help = "BAM file path", metavar = "FILE"),
  make_option(c("--output_dir"), type = "character", help = "Output directory", metavar = "DIR")
)
opt <- parse_args(OptionParser(option_list = option_list))
bam_file <- opt$bam
output_dir <- opt$output_dir
ref <- "/to/path/GRCh38.primary_assembly.genome.fa"
annotation <- readRDS("/to/path/Bambu_transcript_annotations.rds")
se <- bambu(reads = bam_file,annotations = annotation,genome = ref)
writeBambuOutput(se, path = output_dir)
cat("Finished processing:", bam_file, "\n")