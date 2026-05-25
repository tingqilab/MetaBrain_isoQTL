suppressMessages({
  library(coloc)
  library(data.table)
  library(optparse)
})

# ------------------ parameter setting ------------------
option_list <- list(
  make_option("--qtl_gwas", type="character", help="input QTL-GWAS file"),
  make_option("--out_results", type="character", help="coloc all output"),
  make_option("--out_genes", type="character", help="coloc risk gene output"),
  make_option("--pp4", type="numeric", default=0.8, help="PP.H4.abf cutoff")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$qtl_gwas) | is.null(opt$out_results) | is.null(opt$out_genes)) {
  stop("缺少必要参数，请用 --help 查看使用方法")
}

# ------------------ check output file ------------------ #
if (file.exists(opt$out_results)) {
  stop(paste("检测到结果文件已存在:", opt$out_results, "，停止运行以避免覆盖。"))
}

# ------------------ read data ------------------ #
qtl_gwas_df <- fread(opt$qtl_gwas, header = TRUE, data.table = FALSE)
phenotypes <- unique(qtl_gwas_df$phenotype_id)
results_list <- list()

# ------------------ colocalization ------------------ #
for (ph in phenotypes) {
  sub_df <- subset(qtl_gwas_df, phenotype_id == ph)
  # drop NA row
  sub_df$p <- as.numeric(sub_df$P)
  sub_df$pval_nominal <- as.numeric(sub_df$pval_nominal)
  sub_df$b <- as.numeric(sub_df$b)
  sub_df$se <- as.numeric(sub_df$se)
  sub_df$slope <- as.numeric(sub_df$slope)
  sub_df$slope_se <- as.numeric(sub_df$slope_se)
  sub_df$af <- as.numeric(sub_df$af)
  sub_df$af <- as.numeric(sub_df$freq) 
  sub_df <- sub_df[complete.cases(sub_df[, c("p","pval_nominal","b","se","slope","slope_se","af",'freq')]), ]
  sub_df <- sub_df[!duplicated(sub_df$SNP), ] 
  sub_df <- sub_df[sub_df$se > 0 & sub_df$slope_se > 0, ]
  # skipping if no adequet overlapping SNPs
  if (nrow(sub_df) < 10) {
    next
  }
# 假设您的数据是定量性状 (quantitative trait)，所以type设为 "quant"
data1 = list(pvalues = sub_df$p,
             N       = sub_df$N,
             MAF     = sub_df$af,
             beta    = sub_df$b,
             varbeta = sub_df$se^2,
             type    = "quant",          # <--- 将值改为 "quant" 或 "cc"
             snp     = sub_df$SNP)

# data2 保持相似结构，注意 varbeta 的命名和类型
data2 = list(pvalues = sub_df$pval_nominal,
             N       = 2865,             # <--- N 应该是单个值或与 pvalues 等长
             MAF     = sub_df$freq,
             beta    = sub_df$slope,
             varbeta = sub_df$slope_se^2,
             type    = "quant",          # <--- 将值改为 "quant" 或 "cc"
             snp     = sub_df$SNP)

# 注意：coloc要求列表元素必须使用命名方式创建
# 这样做就无需再使用 colnames(data1) = ... 来设置名称了

coloc_res = coloc.abf(data1, data2)
results_list[[ph]] <- coloc_res$summary
}


# ------------------ merge result ------------------
results_df <- do.call(rbind, lapply(names(results_list), function(ph){
  cbind(phenotype_id = ph, t(results_list[[ph]]))
}))
results_df <- as.data.frame(results_df)

# 类型转换
results_df$nsnps <- as.integer(results_df$nsnps)
results_df$PP.H0.abf <- as.numeric(results_df$PP.H0.abf)
results_df$PP.H1.abf <- as.numeric(results_df$PP.H1.abf)
results_df$PP.H2.abf <- as.numeric(results_df$PP.H2.abf)
results_df$PP.H3.abf <- as.numeric(results_df$PP.H3.abf)
results_df$PP.H4.abf <- as.numeric(results_df$PP.H4.abf)

# ------------------ select hits colocalization ------------------
coloc_hits <- subset(results_df, PP.H4.abf > opt$pp4)
coloc_genes <- unique(coloc_hits$phenotype_id)

# ------------------ save result ------------------
write.table(results_df, file=opt$out_results, sep="\t", quote=FALSE, row.names=FALSE)
write.table(coloc_genes, file=opt$out_genes, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

cat("共定位基因数:", length(coloc_genes), "\n")

