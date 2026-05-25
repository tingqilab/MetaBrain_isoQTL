import pandas as pd
import numpy as np
import argparse
import os

def replace_individualID(cohort,df,cov_df,order):
    if cohort=='ROSMAP':
       map_df = cov_df[['individualID','fileID']]
       map_df=map_df[map_df['individualID'].isin(order['IID'])]
       map_df_ordered=map_df.set_index('individualID').loc[order['IID']].reset_index()
       rename_dict = dict(zip(map_df_ordered["fileID"], map_df_ordered["individualID"]))
       df= df.rename(columns=rename_dict)
    if cohort in ['CMC','BrainGVEx','LIBD','HBCC']:
       map_df = cov_df[['individualID','name']]
       map_df=map_df[map_df['individualID'].isin(order['IID'])]
       map_df_ordered=map_df.set_index('individualID').loc[order['IID']].reset_index()
       rename_dict = dict(zip(map_df["name"], map_df["individualID"]))
       df= df.rename(columns=rename_dict)
    if cohort=='Mayo':
       map_df=cov_df[['IID','specimenID']]
       map_df=map_df[map_df['IID'].isin(order['IID'])]
       map_df_ordered=map_df.set_index('IID').loc[order['IID']].reset_index()
       rename_dict=dict(zip(map_df["specimenID"], map_df["IID"]))
       df=df.rename(columns=rename_dict)
    if cohort in ['BM10','BM22','BM36','BM44']:
       map_df=cov_df[['individualIdentifier','sampleIdentifier']]
       map_df=map_df[map_df['individualIdentifier'].isin(order['IID'])]
       map_df_ordered=map_df.set_index('individualIdentifier').loc[order['IID']].reset_index()
       rename_dict=dict(zip(map_df["sampleIdentifier"], map_df["individualIdentifier"]))
       df=df.rename(columns=rename_dict)
    # reorder columns
    ordered_IIDs = order['IID'].tolist()
    ordered_cols = [c for c in ordered_IIDs if c in df.columns]
    other_cols = [c for c in df.columns if c not in ordered_cols]
    return df[other_cols + ordered_cols]


def main(cohort,annotation,feature,measure):
    ## update the sample name to individual ID, and make a bed file by adding chr,start,end.
    transcript_TSS_info='dtranscript_TSS_info.txt'
    Augment_gene_TSS_info='dAugment_gene_TSS_info.txt'
    Gencode_gene_TSS_info='dgene_TSS_info.txt'
    transcript_TSS=pd.read_csv(transcript_TSS_info,sep='\t')
    Gencode_gene_TSS=pd.read_csv(Gencode_gene_TSS_info,sep='\t')
    Augment_gene_TSS=pd.read_csv(Augment_gene_TSS_info,sep='\t')
    Augment_gene_TSS['start']=Augment_gene_TSS['gene_TSS']
    Augment_gene_TSS['end']=Augment_gene_TSS['gene_TSS']+1
    Gencode_gene_TSS['start']=Gencode_gene_TSS['gene_TSS']
    Gencode_gene_TSS['end']=Gencode_gene_TSS['gene_TSS']+1
    transcript_TSS['start']=transcript_TSS['TSS']
    transcript_TSS['end']=transcript_TSS['TSS']+1
    print(f'process cohort {cohort}')
    os.makedirs(f'03_metabrain_cohort_data/{cohort}/phenotype/filtered_result_bed', exist_ok=True)
    feature_id=f'{feature}_id'
    cov_df_file=f'03_metabrain_cohort_data/{cohort}/metadata/{cohort}_cov.txt'
    order_file=f'03_metabrain_cohort_data/{cohort}/covariates/{cohort}_IID_order.txt'
    phenotype_file=f'03_metabrain_cohort_data/{cohort}/phenotype/filtered_result/{cohort}_{annotation}_{feature}_{measure}.txt'
    cov_df=pd.read_csv(cov_df_file,sep='\t')
    order=pd.read_csv(order_file,sep='\t')
    if feature=='transcript':
        tss_info=transcript_TSS[['chrom','transcript_id','start','end']] 
    if feature=='gene':
        tss_info=Augment_gene_TSS[['chrom','gene_id','start','end']]
    tss_info['chrom'] = tss_info['chrom'].str.replace('^chr', '', regex=True)
    if os.path.exists(phenotype_file):
        phenotype=pd.read_csv(phenotype_file,sep='\t')
        phenotype_rename=replace_individualID(cohort,phenotype,cov_df,order)
        phenotype_rename_bed=pd.merge(phenotype_rename,tss_info,how='left',on=feature_id)
        desired_order=['chrom','start','end',feature_id]
        other_cols = [c for c in phenotype_rename_bed.columns if c not in desired_order]
        phenotype_rename_bed_ordered=phenotype_rename_bed[desired_order+other_cols]
        phenotype_rename_bed_ordered = phenotype_rename_bed_ordered[phenotype_rename_bed_ordered['chrom'].isin([str(i) for i in range(1,23)])]
        phenotype_rename_bed_ordered['chrom'] = pd.to_numeric(phenotype_rename_bed_ordered['chrom'], errors='coerce')
        phenotype_rename_bed_ordered['end'] = pd.to_numeric(phenotype_rename_bed_ordered['end'], errors='coerce')
        phenotype_rename_bed_ordered_sort = phenotype_rename_bed_ordered.sort_values(by=['chrom', 'end'])
        meta_cols = phenotype_rename_bed_ordered_sort.columns[:4]
        iid_cols = [iid for iid in order['IID'] if iid in phenotype_rename_bed_ordered_sort.columns]
        phenotype_filtered = phenotype_rename_bed_ordered_sort[list(meta_cols) + iid_cols]
        phenotype_filtered.to_csv(f'03_metabrain_cohort_data/{cohort}/phenotype/filtered_result_bed/{cohort}_{annotation}_{feature}_{measure}.bed',sep='\t',index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert sample names to individual names and make bed file")
    parser.add_argument("--cohort",required=True,help="cohort name(CMC,ROSMAP et al)")
    parser.add_argument("--annotation",required=True,help="Gencode,Metabrain,Metabrain_Augment")
    parser.add_argument("--feature",required=True,help="gene or transcript")
    parser.add_argument("--measure",required=True,help="abundance,count,correct_abundance,correct_count")
    args = parser.parse_args()
    main(cohort=args.cohort,annotation=args.annotation,feature=args.feature,measure=args.measure)

### usage ###
# python script.py --cohort --annotation --feature --measure