### make splicing related phenotype(isoform ratio), use tpm-filtering transcript tpm matrix
import pandas as pd
import os
import argparse

def main(cohort,annotation):
    gene_transcript_map=pd.read_csv('metabrain_Augment_transcript_sturcture_plot.txt',sep='\t')
    gene_transcript_map=gene_transcript_map[['gene_id','transcript_id']].drop_duplicates() 
    work_dir='03_metabrain_cohort_data'
    cohort_dir=f'{work_dir}/{cohort}/phenotype/filtered_result_bed'
    filtered_transcript_bed=f'{cohort_dir}/{cohort}_{annotation}_transcript_abundance.bed'
    TX_TPM_df=pd.read_csv(filtered_transcript_bed,sep='\t')
    TX_TPM_df=pd.merge(TX_TPM_df,gene_transcript_map,on='transcript_id',how='inner')


    gene_tx_counts = TX_TPM_df['gene_id'].value_counts()
    multi_tx_genes = gene_tx_counts[gene_tx_counts > 1].index
    TX_TPM_df_filter = TX_TPM_df[TX_TPM_df['gene_id'].isin(multi_tx_genes)]
    sample_cols = TX_TPM_df_filter.columns[4:-1]  

    def compute_isoform_fraction(group):

        tpm_sum = group[sample_cols].sum(axis=0)  
        return group[sample_cols].div(tpm_sum, axis=1)  
 
    df_TPM_ration = TX_TPM_df_filter.copy()
    df_TPM_ration[sample_cols] = TX_TPM_df.groupby('gene_id', group_keys=False).apply(compute_isoform_fraction)
    df_TPM_ration = df_TPM_ration.drop(columns=['gene_id'])


    missing_ratio = df_TPM_ration[sample_cols].isna().mean(axis=1)
    df_TPM_ration = df_TPM_ration.loc[missing_ratio <= 0.20]
    df_TPM_ration[sample_cols] = df_TPM_ration[sample_cols].apply(lambda x: x.fillna(x.mean()), axis=1)

    out_put_dir=f'{work_dir}/{cohort}/splicing_phenotype/filtered_result_bed'
    os.makedirs(f'{out_put_dir}', exist_ok=True)
    out_path=f'{out_put_dir}/{cohort}_{annotation}_transcript_ratio.bed'
    df_TPM_ration.to_csv(out_path,sep='\t',index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="make isofrom ratio phenotype bed from filtered isoform abundance bed")
    parser.add_argument("--cohort",required=True,help="cohort name(CMC,ROSMAP et al)")
    parser.add_argument("--annotation",required=True,help="Gencode,Metabrain,Metabrain_Augment")
    args = parser.parse_args()
    main(cohort=args.cohort,annotation=args.annotation)


