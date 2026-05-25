library(coloc)
library(dplyr)
library(data.table)
library(optparse)

# -------------------------------
# 1. 命令行参数解析
# -------------------------------
option_list <- list(
  make_option(c("-o", "--output"), type="character", default=NULL,help="Output file path", metavar="character"),
  make_option(c("-p", "--pairs"), type="character", default=NULL,help="Probe pairs file", metavar="character"),
  make_option(c("-i", "--isoform"), type="character", default=NULL, help="sig isoform isoQTLs summary", metavar="character"),
  make_option(c("-e", "--epigenome"), type="character", default=NULL,help="sig epi Probe epiQTLs summary", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

OUT_FILE <- opt$output
pair_path <- opt$pairs
iso_summary_path <- opt$isoform
epi_summary_path <- opt$epigenome


N_iso <-2865
epi_N <- c(
  "caQTL" = 272,
  "brain_mQTL" = 1160,
  "blood_mQTL" = 3701
)

# -------------------------------
# coloc function
# -------------------------------
coloc_function <- function(iso_probe, epi_probe,iso_qtl_summary_df,epi_qtl_summary_df) {

  d1 <- iso_qtl_summary_df[Probe == iso_probe]
  d2 <- epi_qtl_summary_df[Probe == epi_probe]

  if (is.null(d1) | is.null(d2)) return(NULL)
  # merge by SNP
  m <- merge(d1, d2, by = "SNP", suffixes = c(".iso", ".epi"))

  if (nrow(m) < 50) return(NULL)
  # allele harmonization
  same <- m$A1.iso == m$A1.epi & m$A2.iso == m$A2.epi
  flip <- m$A1.iso == m$A2.epi & m$A2.iso == m$A1.epi
  keep <- same | flip
  m <- m[keep, ]
  if (nrow(m) < 50) return(NULL)

  # flip beta & Freq
  m$b.epi[flip[keep]] <- -m$b.epi[flip[keep]]
  m$Freq.epi[flip[keep]] <- 1 - m$Freq.epi[flip[keep]]

  # fill missing epi Freq with iso Freq
  m$Freq.epi[is.na(m$Freq.epi)] <- m$Freq.iso[is.na(m$Freq.epi)]
  if (nrow(m) < 50) return(NULL)
  maf <- pmin(m$Freq.iso, 1 - m$Freq.iso)

  # get qtl_type and epi N
  qtl_type <- unique(m$qtl_type)
  if (length(qtl_type) != 1) return(NULL)
  N_epi <- epi_N[qtl_type]

  # prepare coloc datasets
  dataset1 <- list(
    pvalues = m$p.iso,
    beta = m$b.iso,
    varbeta = m$SE.iso^2,
    snp = m$SNP,
    type = "quant",
    N = N_iso,
    MAF = maf
  )

  dataset2 <- list(
    pvalues = m$p.epi,
    beta = m$b.epi,
    varbeta = m$SE.epi^2,
    snp = m$SNP,
    type = "quant",
    N = N_epi,
    MAF = maf
  )
  coloc_res <- coloc::coloc.abf(dataset1, dataset2)
  data.frame(
    iso_probe = iso_probe,
    epi_probe = epi_probe,
    qtl_type = qtl_type,
    nsnp = nrow(m),
    PP0 = coloc_res$summary["PP.H0.abf"],
    PP1 = coloc_res$summary["PP.H1.abf"],
    PP2 = coloc_res$summary["PP.H2.abf"],
    PP3 = coloc_res$summary["PP.H3.abf"],
    PP4 = coloc_res$summary["PP.H4.abf"]
  )
}

df_pair <- read.table(pair_path, header = TRUE, sep = "\t")
df_pair <- df_pair[ , !names(df_pair) %in% "X"]
epi_qtl_summary_df=fread(epi_summary_path,header = TRUE)
iso_qtl_summary_df=fread(iso_summary_path,header = TRUE)

# -------------------------------
# 4. 循环运行
# -------------------------------
num_pairs <- nrow(df_pair)
results_list <- list()
cat(paste0("开始逐个处理，共计 ", num_pairs, " 对...\n"))
for (i in 1:num_pairs) {
  # 提取当前行的 Probe ID
  current_iso <- df_pair$isoform_probe[i]
  current_epi <- df_pair$epi_probe[i]
  # 运行函数
  res <- coloc_function(current_iso, current_epi,iso_qtl_summary_df, epi_qtl_summary_df)
  # 存储结果
  results_list[[i]] <- res
  # 每处理 50 对打印一次进度
  if (i %% 50 == 0 || i == num_pairs) {
    cat(sprintf("\r进度: [%d / %d] (%.1f%%)", i, num_pairs, (i/num_pairs)*100))
  }
}
cat("\n计算完成！\n")

# -------------------------------
# 5. 合并并输出
# -------------------------------
final_df <- rbindlist(results_list)

if (nrow(final_df) > 0) {
  fwrite(final_df, OUT_FILE, sep = "\t", quote = FALSE)
  cat(paste0("结果已保存至: ", OUT_FILE, "\n"))
} else {
  cat("未产生任何有效结果，请检查数据配对情况。\n")
}