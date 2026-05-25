import pandas as pd
import tensorqtl
import torch
from tensorqtl import pgen, cis, trans, post, genotypeio
import argparse
import os

device = torch.device("cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas: {pd.__version__}")

# --- set command line parameters ---
parser = argparse.ArgumentParser(description="Run eQTL analysis")
parser.add_argument("--cohort", required=True, help="Cohort name, e.g., ROSMAP")
parser.add_argument("--Chresidualdir", required=True, help="this dir have chr1...chr22 phenotype residuals bed file")
parser.add_argument("--chr", required=True, help="Chromosome name, e.g., 1,2,3...")
parser.add_argument("--outdir", required=True, help="input path of output dir")
args = parser.parse_args()

cohort = args.cohort
Chresidualdir=args.Chresidualdir
chr = args.chr
output = args.outdir
# --- define paths ---
work_dir='03_metabrain_cohort_data'
plink_prefix_path = f'{work_dir}/{cohort}/genotype_hg38/{cohort}_chr{chr}'
covariates_file = f'{work_dir}/{cohort}/covariates/{cohort}_genotype_fivePCs_tensorqtl_format.txt'
## how to defien bed file
expression_bed = f'{Chresidualdir}/chr{chr}.bed'
phenotype_name = os.path.basename(Chresidualdir)
output_dir=f'{output}/{phenotype_name}'
os.makedirs(output_dir, exist_ok=True)
prefix = f'{output_dir}/{phenotype_name}'


phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
# --- load genotypes ---
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# --- cis-QTL: permutation / empirical p-values ---
cis_perm_file = prefix + f"_permutation_chr{chr}.csv"
cis_q_file = prefix + f"_permutation_with_q_chr{chr}.csv"
if not os.path.exists(cis_perm_file):
    cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df=covariates_df)
    cis_df.to_csv(cis_perm_file)
else:
    print(f"Permutation cis-QTL file exists: {cis_perm_file}, loading...")
    cis_df = pd.read_csv(cis_perm_file, index_col=0)
if not os.path.exists(cis_q_file):
    post.calculate_qvalues(cis_df, fdr=0.05, qvalue_lambda=0.85)
    cis_df.to_csv(cis_q_file)
else:
    print(f"Q-value file exists: {cis_q_file}, skipping...")