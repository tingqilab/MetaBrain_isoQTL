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
output=args.outdir
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

# --- check if phenotype bed file exists ---
if os.path.exists(expression_bed):
    # --- load phenotypes and covariates ---
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

    # --- load genotypes ---
    pr = genotypeio.PlinkReader(plink_prefix_path)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

    # --- cis-QTL: nominal p-values ---
    cis_nominal_file=f'{prefix}.cis_qtl_pairs.{chr}.parquet'
    if not os.path.exists(cis_nominal_file):
        cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix, covariates_df=covariates_df)
    else:
        print(f"Nominal cis-QTL file exists: {cis_nominal_file}, skipping...")
else:
    print(f"Phenotype bed file does not exist: {expression_bed}, skipping...")




