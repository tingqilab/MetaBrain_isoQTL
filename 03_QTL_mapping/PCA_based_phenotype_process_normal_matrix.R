##  read phenotype bed file and perform normalization before calculating the hidden variables
### for count , we perform TMM normalization and RINT
### for abundance, we perform quantile normalizaiton and RINT

### tpm-filtered phenotype expressed raw count -- TMM normalization (edgeR,correct sequncing size)-- RINT(converting to normal distribution) ###
### TPM, or other ratio normalization values(isoform ratio,junction ratio): TPM -- quantile normalization -- RINT ####

### input two file
## a bed file for phenotype (after filtering,no normalizaiton)
## a covariates file(excluding genotype PC)

### input two optional parameters
## --tier. we use BE methods choose a optimal PC to explain the variance of the matrix, but it totaly is depend on the expression matrix
## without any prior knowledge, so it isn't like to be optimal PC for maximizing QTL discovery power, so tier is necessary
## --thread . tradeoff between running speed and memory relied what computer resource you have 

## examples, here input phenotype bed file, the format is 
## chrom, start, end, gene_id(transcript_id,other phenotype_id), samples1, samples2........
##                                                               expression value..........
## and other input is cov file,the format is 
# individualID    PMI     pH      RIN     ageDeath        diagnosis_Schizophrenia sex_M
# Br827   49.0    6.53    8.3     50.865184       0       1
# Br869   70.5    6.16    6.9     46.862337       1       0
# Br925   29.5    6.48    7.0     66.671229       0       0
# Br942   19.5    6.27    8.4     61.915452       0       1
# Br983   25.0    6.37    9.0     43.382433       1       1
# Br1092  18.0    6.13    9.0     48.945898       0       1
# Br1109  22.0    6.31    8.4     59.221333       0       0
# Br1113  51.5    6.53    8.7     51.905596       0       0
# Br1183  29.5    6.72    8.6     26.065053       1       0
# Br1189  23.5    6.49    8.1     70.879421       1       0
# Br1203  32.0    6.62    8.2     24.334684       1       1
# Br1225  22.5    6.38    8.5     32.674405       0       1
# Br1304  64.0    6.87    8.1     58.566969       0       1
# Br1339  89.0    6.65    7.4     60.075566       1       0

# note the sample ID have to matching one by one

## usage generally 
# Rscript PCA_based_phenotype_process_normal_matrix.R \
#   --bed /to/path/test.bed \
#   --cov /to/path/test.txt \
#   --tier=-15,-10,-5,0,5,10,15,20 \   ## there we ues BE methods to define a PC number and get the different tier around the number
#   --threads 8 \
#   --normalization 'raw_count_normas' \
#   --output /to/path/output_dir/
#   --prefix output_prefix


### the memory usage and running time
### we test in a expression matrix 800 samples and 140000 gene , use 4 core, 4 threads , 120G memory , need 25 mins to run
### we test in 8000 samples and 140000 gene , use 4 parallel and  480 threads, 400 G memory , need 80 mins to run



## frist install the https://github.com/heatherjzhou/PCAForQTL
library(limma)  
library(edgeR)
library(DESeq2)
library(glue)
library(optparse)
library(BiocParallel)


option_list <- list(
  make_option(c("--bed"), type="character"),
  make_option(c("--cov"), type="character"),
  make_option(c("--tier"), type="character", default="0"),  # default 0
  make_option(c("--threads"), type = "integer", default=1),  # default 1
  make_option(c("--normalization"), type = "character"), ### two available option, ratio_norma(for ration phenotype,such psi, iso ratio),raw_count_norma
  make_option(c("--output"), type = "character"), ## output_dir
  make_option(c("--prefix"), type = "character") ## output_prefix
)

opt <- parse_args(OptionParser(option_list=option_list))
cat("bed_file =", opt$bed, "\n")
cat("covariates_file =", opt$cov, "\n")
cat("tier_input =", opt$tier, "\n")
cat("threads =", opt$threads, "\n")


## covert string to numeric vector
tier <- as.numeric(strsplit(opt$tier, ",")[[1]])
cat("tier vector =", tier, "\n")

register(MulticoreParam(workers = opt$threads))

Sys.setenv(
  OMP_NUM_THREADS = opt$threads,
  MKL_NUM_THREADS = opt$threads,
  OPENBLAS_NUM_THREADS = opt$threads
)

RINT <- function(x){
  r <- rank(x, na.last = "keep", ties.method = "random")
  qnorm((r - 3/8) / (sum(!is.na(x)) + 1/4))
}

#### quantile normalization + RINT
qn_rint_matrix <- function(expr_mat){
  # expr_mat: row = gene/transcript, col = samples
  expr_qn <- normalizeBetweenArrays(expr_mat, method = "quantile")
  expr_rint <- t(apply(expr_qn, 1, RINT))
  rownames(expr_rint) <- rownames(expr_mat)
  colnames(expr_rint) <- colnames(expr_mat)
  return(expr_rint)
}

ratio_norma <- function(bed){
  bed_info <- bed[, 1:4]              # metadata chrom, start, end, gene_id
  expr_mat <- as.matrix(bed[, -c(1:4)])  # samples 
  rownames(expr_mat) <- bed$gene_id
  # quantile normalization + RINT
  expr_rint <- qn_rint_matrix(expr_mat)
  bed <- cbind(bed_info, expr_rint)
  # ordering
  bed <- bed[order(as.numeric(bed$chrom), bed$start), ]
  return(bed)
}

raw_count_norma <- function(bed){
  # gene info
  gene_info <- bed[, 1:4]
  # raw count matrix
  counts <- as.matrix(bed[, -c(1:4)])
  # edgeR TMM
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  # logCPM
  logcpm <- cpm(dge, log = TRUE, prior.count = 1)
  # RINT per gene
  expr_rint <- t(apply(logcpm, 1, RINT))
  # combine
  bed_rint <- cbind(gene_info, expr_rint)
  return(bed_rint)
}


getResiduals<-function(dataResponse,dataPredictors){
  vecOfOnes<-rep(1,nrow(dataResponse)) 
  X<-cbind(vecOfOnes,as.matrix(dataPredictors)) 
  betaEstimate<-solve(t(X)%*%X)%*%t(X)%*%as.matrix(dataResponse) 
  toReturn<-as.matrix(dataResponse)-X%*%betaEstimate # Residuals.
  return(toReturn)
}


### this function input normalized bed file, covariates excluding genotype PC, optional tier(default is 0)
## alternative tier is a list, such as (-20,-15,-10,-5,0,5,10)
## usually BE is upper selection for PC
PCA_getResiduals<-function(normalized_bed,cov,tier){
    ### 0. PCA process for gene 
    expr<-t(normalized_bed[,-(1:4)])
    metadata<-normalized_bed[,(1:4)]
    prcompResult<-prcomp(expr,center=TRUE,scale.=TRUE) #This should take less than a minute.
    PCs<-prcompResult$x  # The columns are PCs.
    # 1. Elbow plot
    resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult)
    # 2. BE methods, usually it will generate the higher PC number
    RNGkind("L'Ecuyer-CMRG")
    set.seed(1)
    resultRunBE<-PCAForQTL::runBE(expr,B=20,alpha=0.05,mc.cores=opt$threads)
    K_BE<-resultRunBE$numOfPCsChosen
    cat("K_pc =:", K_BE, "\n")
    ## choose K=K_BE for PC
    if(length(tier) == 1 && tier == 0) {
        PCsTop<-PCs[,1:K_BE]
        knownCovariatesFiltered<-PCAForQTL::filterKnownCovariates(cov,PCsTop,unadjustedR2_cutoff=0.9)
        PCsTop<-scale(PCsTop) #Optional. Could be helpful for avoiding numerical inaccuracies.
        covariatesToUse<-cbind(knownCovariatesFiltered,PCsTop)
        ## get Residuals
        Residuals=getResiduals(expr,covariatesToUse)
        final_bed<- cbind(metadata, t(Residuals))
        return(setNames(list(final_bed), paste0("tier_PC", K_BE))) 
    }
    
    # 6. tier >0 多输出
    residual_list <- list()
    for(t in tier){
        nPC <- K_BE + t
        nPC <- max(1, min(nPC, ncol(PCs)))  # 确保在合理范围
        PCsTop <- scale(PCs[,1:nPC])
        knownCovariatesFiltered<-PCAForQTL::filterKnownCovariates(cov,PCsTop,unadjustedR2_cutoff=0.9)
        covariatesToUse <- cbind(knownCovariatesFiltered, PCsTop)
        Residuals <- getResiduals(expr, covariatesToUse)
        residual_list[[paste0("tier_", nPC)]] <- cbind(metadata, t(Residuals))
    }
    return(residual_list)
}

saveResiduals <- function(residual_list,output_dir,prefix){
  for(tier_name in names(residual_list)){
      # tier_name 是 list 的名字，比如 tier_10
      nPC <- sub("tier_", "", tier_name)  # 取出数字
      outfile <- glue("{output_dir}/{prefix}_tier_PC{nPC}.bed")
      write.table(residual_list[[tier_name]], file=outfile, sep="\t", row.names=FALSE, quote=FALSE)
  }
}

## read cov and bed file
cov_path <- opt$cov
cov <- read.delim(cov_path, header = TRUE, row.names = 1, sep = "\t") 
bet_path <- opt$bed
bed_df=read.delim(bet_path,sep='\t',check.names = FALSE)

### normalization
normalization_method <- opt$normalization
if(normalization_method=='ratio_norma'){
    normalize_bed_df=ratio_norma(bed_df)
}

if(normalization_method=='raw_count_norma'){
    normalize_bed_df=raw_count_norma(bed_df)
}

# ---- output dir----
output_dir <- opt$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
### PCA get the Residuals and save the tier result ###
# output dir
output_dir <- opt$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
residual_gene_count_list <- PCA_getResiduals(normalize_bed_df,cov,tier)
saveResiduals(residual_gene_count_list,output_dir,opt$prefix)