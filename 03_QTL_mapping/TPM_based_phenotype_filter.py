### preparing phenotype, filter, normalization, RANT
### define some function
### 1.including transcript which tmp>tmp_cutoff in at least prop_cutoff samples
import os
import pandas as pd
import numpy as np
import argparse

def filter_phenotype(tmp_matrix,TPM_cutoff=0.1, prop_cutoff=0.2):
    prop = (tmp_matrix.iloc[:, 1:] > TPM_cutoff).mean(axis=1)
    filtered = tmp_matrix[prop >= prop_cutoff].copy()
    return filtered
### 2. converting counts to tpm
def calculate_tpm_from_count_length(count_df,length_df):
    transcript_ids = count_df['transcript_id']
    counts = count_df.iloc[:, 1:]
    lengths = length_df.iloc[:, 1:]
    lengths_kb = lengths / 1000
    rpk = counts / lengths_kb
    scaling_factor = rpk.sum(axis=0) / 1e6  # 每列单独算
    tpm = rpk / scaling_factor
    tpm_df = pd.concat([transcript_ids, tpm], axis=1)
    return tpm_df
### 3. rename the sample name to individual ID
def replace_individualID(cohort,df,cov_df,order):
    if cohort=='ROSMAP':
       map_df = cov_df[['individualID','fileID']]
       map_df=map_df[map_df['individualID'].isin(order['IID'])]
       map_df_ordered=map_df.set_index('individualID').loc[order['IID']].reset_index()
       rename_dict = dict(zip(map_df_ordered["fileID"], map_df_ordered["individualID"]))
       df= df.rename(columns=rename_dict)
    if cohort in ['CMC','BrainGVEX','LIBD','HBCC']:
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
    return df[ordered_cols + other_cols]

### stats filter result
def transcript_info_stats(transcript_info,Gencode_transcript_abundance_filter):
    transcript_info_cohort=transcript_info[transcript_info['transcript_id'].isin(Gencode_transcript_abundance_filter['transcript_id'])]
    transcript_info_cohort['novel_gene']=transcript_info_cohort['gene_id'].str.startswith('novel')
    annotated_gene_num=transcript_info_cohort[~transcript_info_cohort['novel_gene']]['gene_id'].nunique()
    novel_gene_num=transcript_info_cohort[transcript_info_cohort['novel_gene']]['gene_id'].nunique()
    annotated_transcript_num=len(transcript_info_cohort[~transcript_info_cohort['novel_transcript']])
    novel_transcript_num=len(transcript_info_cohort[transcript_info_cohort['novel_transcript']])
    return transcript_info_cohort,annotated_gene_num,novel_gene_num,annotated_transcript_num,novel_transcript_num

### filtering 
def find_best_cutoff_for_Augment_transcript(goal_annotated_transcript_num,Metabrain_Augment_transcript_abundance):
    cutoffs = np.arange(0.20, 0.50 + 0.01, 0.01)
    records = []
    for cutoff in cutoffs:
        transcript_abundance_filter = filter_phenotype(Metabrain_Augment_transcript_abundance,TPM_cutoff=0.1,prop_cutoff=cutoff)
        annotated_transcript_num = len(transcript_abundance_filter[transcript_abundance_filter['transcript_id'].str.startswith('ENST')])
        records.append({'cutoff': cutoff,'Augment_annotated_transcript_num': annotated_transcript_num,'diff': abs(annotated_transcript_num - goal_annotated_transcript_num)})
    result_df = pd.DataFrame(records)
    # find cutoff
    best_row = result_df.loc[result_df['diff'].idxmin()]
    best_Augment_transcript = filter_phenotype(Metabrain_Augment_transcript_abundance, 0.1, best_row['cutoff'])
    print(f"cutoff for Augment transcript is {best_row['cutoff']}, diff is {best_row['diff']}")
    return best_Augment_transcript

def find_best_cutoff_for_Gencode_gene(goal_annotated_gene_num,Metabrain_Gencode_gene_abundance):
    cutoffs = np.arange(0.20, 0.50 + 0.01, 0.01)
    records = []
    for cutoff in cutoffs:
        gene_abundance_filter = filter_phenotype(Metabrain_Gencode_gene_abundance,TPM_cutoff=0.1,prop_cutoff=cutoff)
        annotated_gene_num = len(gene_abundance_filter)
        records.append({'cutoff': cutoff,'annotated_gene_num': annotated_gene_num,'diff': abs(annotated_gene_num - goal_annotated_gene_num)})
    result_df = pd.DataFrame(records)
    # find cutoff
    best_row = result_df.loc[result_df['diff'].idxmin()]
    best_Gencode_gene = filter_phenotype(Metabrain_Gencode_gene_abundance, 0.1, best_row['cutoff'])
    print(f"cutoff for Gencode gene is {best_row['cutoff']}, diff is {best_row['diff']}")
    return best_Gencode_gene


def main(cohort):
    ### 1.TPM-based-filtering
    ## define path per cohorts
    work_dir='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data'
    os.makedirs(f'{work_dir}/{cohort}/phenotype/filtered_result',exist_ok=True)
    out_put_dir=f'{work_dir}/{cohort}/phenotype/filtered_result'
    cohort_work_dir=f'{work_dir}/{cohort}/phenotype/merge_result/{cohort}'
    cov_df_file=f'{work_dir}/{cohort}/metadata/{cohort}_cov.txt'
    order_file=f'{work_dir}/{cohort}/covariates/{cohort}_IID_order.txt'

    Gencode_gene_abundance_path=f'{cohort_work_dir}_Gencode_gene_abundance.txt'
    Gencode_gene_counts_path=f'{cohort_work_dir}_Gencode_gene_counts.txt'
    Gencode_transcript_abundance_path=f'{cohort_work_dir}_Gencode_transcript_abundance.txt'
    Gencode_transcript_counts_path=f'{cohort_work_dir}_Gencode_transcript_counts.txt'
    Gencode_transcript_length_path=f'{cohort_work_dir}_Gencode_transcript_length.txt'
    Gencode_transcript_correct_counts_path=f'{cohort_work_dir}_Gencode_transcript_RTAcorrect_count.txt'

    Metabrain_Augment_gene_abundance_path=f'{cohort_work_dir}_Metabrain_Augment_gene_abundance.txt'
    Metabrain_Augment_gene_counts_path=f'{cohort_work_dir}_Metabrain_Augment_gene_counts.txt'
    Metabrain_Augment_transcript_abundance_path=f'{cohort_work_dir}_Metabrain_Augment_transcript_abundance.txt'
    Metabrain_Augment_transcript_counts_path=f'{cohort_work_dir}_Metabrain_Augment_transcript_counts.txt'
    Metabrain_Augment_transcript_length_path=f'{cohort_work_dir}_Metabrain_Augment_transcript_length.txt'
    Metabrain_Augment_transcript_correct_counts_path=f'{cohort_work_dir}_Metabrain_Augment_transcript_RTAcorrect_count.txt'

    Metabrain_gene_abundance_path=f'{cohort_work_dir}_Metabrain_gene_abundance.txt'
    Metabrain_gene_counts_path=f'{cohort_work_dir}_Metabrain_gene_counts.txt'
    Metabrain_transcript_abundance_path=f'{cohort_work_dir}_Metabrain_transcript_abundance.txt'
    Metabrain_transcript_counts_path=f'{cohort_work_dir}_Metabrain_transcript_counts.txt'
    Metabrain_transcript_length_path=f'{cohort_work_dir}_Metabrain_transcript_length.txt'
    Metabrain_transcript_correct_counts_path=f'{cohort_work_dir}_Metabrain_transcript_RTAcorrect_count.txt'

    transcript_TSS_info='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/transcript_TSS_info/transcript_TSS_info.txt'
    Augment_gene_TSS_info='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/transcript_TSS_info/Augment_gene_TSS_info.txt'
    Gencode_gene_TSS_info='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/transcript_TSS_info/gene_TSS_info.txt'
    transcript_info=pd.read_csv('/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/01_construct_transcript_panel/transcript_final_panel_analysis/metabrain_Augment_transcript_sturcture_plot.txt',sep='\t')
    transcript_info=transcript_info[transcript_info['feature_type']=='transcript']

    ### 1.process correct transcript count to abundance
    Gencode_transcript_correct_counts=pd.read_csv(Gencode_transcript_correct_counts_path,sep='\t')
    Gencode_transcript_length=pd.read_csv(Gencode_transcript_length_path,sep='\t')
    Gencode_transcript_correct_abundance=calculate_tpm_from_count_length(Gencode_transcript_correct_counts,Gencode_transcript_length)
    Gencode_transcript_correct_abundance.to_csv(f'{cohort_work_dir}_Gencode_transcript_RTAcorrect_abundance.txt',sep='\t',index=False)
    Metabrain_transcript_correct_counts=pd.read_csv(Metabrain_transcript_correct_counts_path,sep='\t')
    Metabrain_transcript_length=pd.read_csv(Metabrain_transcript_length_path,sep='\t')
    Metabrain_transcript_correct_abundance=calculate_tpm_from_count_length(Metabrain_transcript_correct_counts,Metabrain_transcript_length)
    Metabrain_transcript_correct_abundance.to_csv(f'{cohort_work_dir}_Metabrain_transcript_RTAcorrect_abundance.txt',sep='\t',index=False)
    Metabrain_Augment_transcript_correct_counts=pd.read_csv(Metabrain_Augment_transcript_correct_counts_path,sep='\t')
    Metabrain_Augment_transcript_length=pd.read_csv(Metabrain_Augment_transcript_length_path,sep='\t')
    Metabrain_Augment_transcript_correct_abundance=calculate_tpm_from_count_length(Metabrain_Augment_transcript_correct_counts,Metabrain_Augment_transcript_length)
    Metabrain_Augment_transcript_correct_abundance.to_csv(f'{cohort_work_dir}_Metabrain_Augment_transcript_RTAcorrect_abundance.txt',sep='\t',index=False)

    ## 2. TPM-based filtering, 
    # transcript_filter: fix tpm>0.1 in at least half of samples in Gencode, then adjust the more release cutoff in Augment annotation
    ### for panel Metabrain, here we give a assumption , that is the transcripts belong to Metabrain is brain expressed transcript, or we can call metabrain is
    # brain specific full length transcript panel
    # gene_filter: fix tpm>0.1 in at least 20% samples in Augment , then adjust gencode 

    Gencode_transcript_abundance=pd.read_csv(Gencode_transcript_abundance_path,sep='\t')
    Metabrain_transcript_abundance=pd.read_csv(Metabrain_transcript_abundance_path,sep='\t')
    Metabrain_Augment_transcript_abundance=pd.read_csv(Metabrain_Augment_transcript_abundance_path,sep='\t')
    Gencode_transcript_abundance_filter=filter_phenotype(Gencode_transcript_abundance,TPM_cutoff=0.1,prop_cutoff=0.5)
    Metabrain_transcript_abundance_filter=filter_phenotype(Metabrain_transcript_abundance,TPM_cutoff=0.1,prop_cutoff=0.5)
    goal_annotated_transcript_num=len(Gencode_transcript_abundance_filter)

    Gencode_gene_abundance=pd.read_csv(Gencode_gene_abundance_path,sep='\t')
    Metabrain_gene_abundance=pd.read_csv(Metabrain_gene_abundance_path,sep='\t')
    Metabrain_Augment_gene_abundance=pd.read_csv(Metabrain_Augment_gene_abundance_path,sep='\t')
    Metabrain_gene_abundance_filter=filter_phenotype(Metabrain_gene_abundance,TPM_cutoff=0.1,prop_cutoff=0.2)
    Metabrain_Augment_gene_abundance_filter=filter_phenotype(Metabrain_Augment_gene_abundance,TPM_cutoff=0.1,prop_cutoff=0.2)
    goal_annotated_gene_num=len(Metabrain_Augment_gene_abundance_filter[(Metabrain_Augment_gene_abundance_filter['gene_id'].str.startswith('EN'))&(Metabrain_Augment_gene_abundance_filter['gene_id'].str.len()<25)])

    ### adjust cutoff 
    Gencode_gene_abundance_filter=find_best_cutoff_for_Gencode_gene(goal_annotated_gene_num,Gencode_gene_abundance) #prop_cutoff>0.2
    Metabrain_Augment_transcript_abundance_filter=find_best_cutoff_for_Augment_transcript(goal_annotated_transcript_num,Metabrain_Augment_transcript_abundance) #prop_cutoff < 0.5

    ### get gene count and transcript correct abundance for downstrem analysis
    # Gencode_transcript_correct_abundance_path=f'{cohort_work_dir}_Gencode_transcript_RTAcorrect_abundance.txt'
    # Metabrain_transcript_correct_abundance_path=f'{cohort_work_dir}_Metabrain_transcript_RTAcorrect_abundance.txt'
    # Metabrain_Augment_transcript_correct_abundance_path=f'{cohort_work_dir}_Metabrain_Augment_transcript_RTAcorrect_abundance.txt'
    # Gencode_transcript_correct_abundance=pd.read_csv(Gencode_transcript_correct_abundance_path,sep='\t')
    # Metabrain_transcript_correct_abundance=pd.read_csv(Metabrain_transcript_correct_abundance_path,sep='\t')
    # Metabrain_Augment_transcript_correct_abundance=pd.read_csv(Metabrain_Augment_transcript_correct_abundance_path,sep='\t')

    Gencode_gene_counts=pd.read_csv(Gencode_gene_counts_path,sep='\t')
    Metabrain_gene_counts=pd.read_csv(Metabrain_gene_counts_path,sep='\t')
    Metabrain_Augment_gene_counts=pd.read_csv(Metabrain_Augment_gene_counts_path,sep='\t')

    ### save the filtered file to dir
    Gencode_transcript_correct_abundance_filter=Gencode_transcript_correct_abundance[Gencode_transcript_correct_abundance['transcript_id'].isin(Gencode_transcript_abundance_filter['transcript_id'])]
    Metabrain_transcript_correct_abundance_filter=Metabrain_transcript_correct_abundance[Metabrain_transcript_correct_abundance['transcript_id'].isin(Metabrain_transcript_abundance_filter['transcript_id'])]
    Metabrain_Augment_transcript_correct_abundance_filter=Metabrain_Augment_transcript_correct_abundance[Metabrain_Augment_transcript_correct_abundance['transcript_id'].isin(Metabrain_Augment_transcript_abundance_filter['transcript_id'])]
    Gencode_gene_counts_filter=Gencode_gene_counts[Gencode_gene_counts['gene_id'].isin(Gencode_gene_abundance_filter['gene_id'])]
    Metabrain_gene_counts_filter=Metabrain_gene_counts[Metabrain_gene_counts['gene_id'].isin(Metabrain_gene_abundance_filter['gene_id'])]
    Metabrain_Augment_gene_counts_filter=Metabrain_Augment_gene_counts[Metabrain_Augment_gene_counts['gene_id'].isin(Metabrain_Augment_gene_abundance_filter['gene_id'])]
    
    Gencode_transcript_counts=pd.read_csv(Gencode_transcript_counts_path,sep='\t')
    Gencode_transcript_counts_filter=Gencode_transcript_counts[Gencode_transcript_counts['transcript_id'].isin(Gencode_transcript_abundance_filter['transcript_id'])]
    Gencode_transcript_correct_counts_filter= Gencode_transcript_correct_counts[Gencode_transcript_correct_counts['transcript_id'].isin(Gencode_transcript_abundance_filter['transcript_id'])]
    
    Metabrain_transcript_counts=pd.read_csv(Metabrain_transcript_counts_path,sep='\t')
    Metabrain_transcript_counts_filter=Metabrain_transcript_counts[Metabrain_transcript_counts['transcript_id'].isin(Metabrain_transcript_abundance_filter['transcript_id'])]
    Metabrain_transcript_correct_counts_filter= Metabrain_transcript_correct_counts[Metabrain_transcript_correct_counts['transcript_id'].isin(Metabrain_transcript_abundance_filter['transcript_id'])]
    
    Metabrain_Augment_transcript_counts=pd.read_csv(Metabrain_Augment_transcript_counts_path,sep='\t')
    Metabrain_Augment_transcript_counts_filter=Metabrain_Augment_transcript_counts[Metabrain_Augment_transcript_counts['transcript_id'].isin(Metabrain_Augment_transcript_abundance_filter['transcript_id'])]
    Metabrain_Augment_transcript_correct_counts_filter= Metabrain_Augment_transcript_correct_counts[Metabrain_Augment_transcript_correct_counts['transcript_id'].isin(Metabrain_Augment_transcript_abundance_filter['transcript_id'])] 


    Gencode_transcript_abundance_filter.to_csv(f'{out_put_dir}/{cohort}_Gencode_transcript_abundance.txt',sep='\t',index=False)
    Gencode_transcript_correct_abundance_filter.to_csv(f'{out_put_dir}/{cohort}_Gencode_transcript_correct_abundance.txt',sep='\t',index=False)
    Gencode_gene_abundance_filter.to_csv(f'{out_put_dir}/{cohort}_Gencode_gene_abundance.txt',sep='\t',index=False)
    Gencode_gene_counts_filter.to_csv(f'{out_put_dir}/{cohort}_Gencode_gene_count.txt',sep='\t',index=False)
    Gencode_transcript_counts_filter.to_csv(f'{out_put_dir}/{cohort}_Gencode_transcript_count.txt',sep='\t',index=False)
    Gencode_transcript_correct_counts_filter.to_csv(f'{out_put_dir}/{cohort}_Gencode_transcript_correct_count.txt',sep='\t',index=False)

    Metabrain_transcript_abundance_filter.to_csv(f'{out_put_dir}/{cohort}_Metabrain_transcript_abundance.txt',sep='\t',index=False)
    Metabrain_transcript_correct_abundance_filter.to_csv(f'{out_put_dir}/{cohort}_Metabrain_transcript_correct_abundance.txt',sep='\t',index=False)
    Metabrain_gene_abundance_filter.to_csv(f'{out_put_dir}/{cohort}_Metabrain_gene_abundance.txt',sep='\t',index=False)
    Metabrain_gene_counts_filter.to_csv(f'{out_put_dir}/{cohort}_Metabrain_gene_count.txt',sep='\t',index=False)
    Metabrain_transcript_counts_filter.to_csv(f'{out_put_dir}/{cohort}_Metabrain_transcript_count.txt',sep='\t',index=False)
    Metabrain_transcript_correct_counts_filter.to_csv(f'{out_put_dir}/{cohort}_Metabrain_transcript_correct_count.txt',sep='\t',index=False)
    
    Metabrain_Augment_transcript_abundance_filter.to_csv(f'{out_put_dir}/{cohort}_Metabrain_Augment_transcript_abundance.txt',sep='\t',index=False)
    Metabrain_Augment_transcript_correct_abundance_filter.to_csv(f'{out_put_dir}/{cohort}_Metabrain_Augment_transcript_correct_abundance.txt',sep='\t',index=False)
    Metabrain_Augment_gene_abundance_filter.to_csv(f'{out_put_dir}/{cohort}_Metabrain_Augment_gene_abundance.txt',sep='\t',index=False)
    Metabrain_Augment_gene_counts_filter.to_csv(f'{out_put_dir}/{cohort}_Metabrain_Augment_gene_count.txt',sep='\t',index=False)
    Metabrain_Augment_transcript_counts_filter.to_csv(f'{out_put_dir}/{cohort}_Metabrain_Augment_transcript_count.txt',sep='\t',index=False)
    Metabrain_Augment_transcript_correct_counts_filter.to_csv(f'{out_put_dir}/{cohort}_Metabrain_Augment_transcript_correct_count.txt',sep='\t',index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="filter phenotype")
    parser.add_argument("--cohort",required=True,help="cohort name")
    args = parser.parse_args()
    main(cohort=args.cohort)

## usage
# python your_script.py --cohort ROSMAP