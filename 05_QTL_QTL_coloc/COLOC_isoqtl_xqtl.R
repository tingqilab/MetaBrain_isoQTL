suppressMessages({
  library(coloc)
  library(data.table)
  library(optparse)
})

# ------------------ parameter setting ------------------
option_list <- list(
  make_option("--qtl_qtl", type="character", help="input QTL isoform summary file for one gene"),
  make_option("--out_results", type="character", help="coloc all output"),
  make_option("--out_genes", type="character", help="coloc hit isoform pairs"),
  make_option("--pp4", type="numeric", default=0.8, help="PP.H4.abf cutoff")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$qtl_qtl) | is.null(opt$out_results) | is.null(opt$out_genes)) {
  stop("缺少必要参数，请用 --help 查看使用方法")
}

# ------------------ check output file ------------------ #
if (file.exists(opt$out_results)) {
  stop(paste("检测到结果文件已存在:", opt$out_results, "，停止运行以避免覆盖。"))
}
# ------------------ read data ------------------ #
qtl_df <- fread(opt$qtl_qtl, header = TRUE, data.table = FALSE)


# ------------------ isoform list ------------------ #
Probes <- unique(qtl_df$Probe)      
isoforms <- Probes[grepl("^(ENST|TCONS)", Probes)]  ### isoform expression Probe 
non_isoforms <- Probes[!grepl("^(ENST|TCONS)", Probes)]  ### non isoform expression Probe
gene_name <- unique(qtl_df$Gene)      ### 可用来标注最终输出
results_list <- list()

# =================================================================
#             对同一基因的 isoform 两两共定位
# =================================================================

# 获取两两所有 isoform 组合
pairs <- expand.grid(isoforms, non_isoforms, stringsAsFactors = FALSE)
# print(pairs)
for (i in seq_len(nrow(pairs))) {
  iso1 <- pairs[i, 1]
  iso2 <- pairs[i, 2]  ### non isoform expression Probe
  print(c(iso1, iso2))
  df1 <- subset(qtl_df, Probe == iso1)
  df2 <- subset(qtl_df, Probe == iso2)
  # 只保留两个 isoform 共享的 SNP
  merged <- merge(df1, df2, by="SNP", suffixes = c("_1", "_2"))
  merged <- subset(merged, complete.cases(merged[, c("p_1","p_2","Freq_1","Freq_2")]))

  if (nrow(merged) < 10) next

  # coloc 输入
  data1 <- list(
    pvalues = merged$p_1,
    N = 2865,
    MAF = merged$Freq_1,
    beta = merged$b_1,
    varbeta = merged$SE_1^2,
    type = "quant",
    snp = merged$SNP
  )
  data2 <- list(
    pvalues = merged$p_2,
    N = 2865,
    MAF = merged$Freq_2,
    beta = merged$b_2,
    varbeta = merged$SE_2^2,
    type = "quant",
    snp = merged$SNP
  )

  stopifnot(identical(data1$snp, data2$snp))
  coloc_res <- coloc.abf(data1, data2)
  ph <- paste(iso1, iso2, sep = "_VS_")    ### <-- 修改，输出为 iso1_VS_iso2
  results_list[[ph]] <- coloc_res$summary
}

# ------------------ merge result ------------------
results_df <- rbindlist(lapply(names(results_list), function(ph){
  cbind(iso_pair = ph, gene = gene_name, as.data.frame(t(results_list[[ph]])))
}), fill = TRUE)

# 类型转换
cols <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
results_df[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
# ------------------ select hits ------------------
coloc_hits <- subset(results_df, PP.H4.abf > opt$pp4)
# ------------------ save result ------------------
dir.create(dirname(opt$out_results), showWarnings = FALSE, recursive = TRUE)
fwrite(results_df, opt$out_results, sep="\t", quote=FALSE)
fwrite(coloc_hits, opt$out_genes, sep="\t", quote=FALSE)
cat("共定位 Probe pairs 数:", nrow(coloc_hits), "\n")