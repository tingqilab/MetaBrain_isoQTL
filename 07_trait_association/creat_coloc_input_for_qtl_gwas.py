import pandas as pd
import os
import argparse


def merge_gwas_qtl(qtl_df, gwas_df):
    gwas_qtl = pd.merge(
        gwas_df,
        qtl_df,
        how='inner',
        on='SNP',
        suffixes=('_gwas', '_qtl')
    )
    return gwas_qtl


def transform_freq(gwas_eqtl):
    mask = (
        (gwas_eqtl['A1_gwas'] == gwas_eqtl['A2_qtl']) &
        (gwas_eqtl['A2_gwas'] == gwas_eqtl['A1_qtl'])
    )

    gwas_eqtl.loc[mask, 'b'] = -gwas_eqtl.loc[mask, 'b']
    gwas_eqtl.loc[mask, ['A1_gwas', 'A2_gwas']] = \
        gwas_eqtl.loc[mask, ['A2_gwas', 'A1_gwas']].values
    gwas_eqtl.loc[mask, 'freq'] = 1 - gwas_eqtl.loc[mask, 'freq']

    return gwas_eqtl


def process_gwas_qtl_pairs(gwas_df, qtl_df, output_path):
    gwas_qtl = merge_gwas_qtl(qtl_df, gwas_df)
    gwas_qtl_summary = transform_freq(gwas_qtl)

    gwas_qtl_summary_clean = gwas_qtl_summary[
        (gwas_qtl_summary['freq'] - gwas_qtl_summary['af']).abs() < 0.1
    ]

    print(
        f"{gwas_qtl_summary_clean['SNP'].nunique()}/"
        f"{gwas_qtl_summary['SNP'].nunique()} SNPs with consistent allele frequency"
    )

    gwas_qtl_summary_clean.to_csv(output_path, sep='\t', index=False)

    return gwas_qtl_summary_clean


def main():
    parser = argparse.ArgumentParser(
        description="Merge GWAS locus file with QTL summary for coloc input"
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory"
    )

    parser.add_argument(
        "--gwas_locus",
        required=True,
        help="GWAS locus file"
    )

    parser.add_argument(
        "--qtl_summary",
        required=True,
        help="QTL summary file"
    )

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    gwas_name = args.gwas_locus.split('/')[-2]
    locus_label = args.gwas_locus.split('/')[-1].split('.')[0]
    qtl_name = args.qtl_summary.split('/')[-1].split('.')[0]

    gwas_locus_df = pd.read_csv(
        args.gwas_locus,
        sep='\t',
        header=1,
        names=['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'P', 'N', 'chr', 'pos']
    )

    qtl_summary_df = pd.read_csv(
        args.qtl_summary,
        sep='\t',
        usecols=[0, 1, 2, 3, 4, 5, 6, 9, 11, 12, 13],
        header=1,
        names=[
            'SNP', 'chr', 'pos', 'A1', 'A2', 'af',
            'phenotype_id', 'gene',
            'slope', 'slope_se', 'pval_nominal'
        ]
    )
    output_path = f"{args.output_dir}/{gwas_name}-{locus_label}-{qtl_name}-coloc_input.txt"
    print('result save to ', output_path)
    process_gwas_qtl_pairs(
        gwas_locus_df,
        qtl_summary_df,
        output_path
    )

if __name__ == "__main__":
    main()


## usage example:
# python merge_gwas_qtl.py \
#   --output_dir /qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/qtl_clumping_gwas_coloc_input \
#   --gwas_locus /qitinglab/sharedata/GWAS_summary/brain_disease_gwas_newest_EUR/Clumping_GWAS_QtlV_split_chr_for_coloc/AD_2022_NatGenet/locus_chr4_460186_1460186.csv \
#   --qtl_summary /qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_all/Gencode_isoQTL_chr4.txt