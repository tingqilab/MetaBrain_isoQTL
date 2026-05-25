import pandas as pd
import os
import sys
from collections import defaultdict
from itertools import product
chr = sys.argv[1]
out_dir='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/isoQTL_xQTL_input'
qtl_summary_dir = '/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_all'
qtl_lead_dir = '/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_per_probe_lead'
epiQTL_sigProbe_summary_data='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_unique_isoQTL_genetic_control/epi_sigQTL_data'
epiQTltypes=['caQTL','brain_mQTL','blood_mQTL']

## read isoQTl
isoQTL_all=pd.read_csv(f'{qtl_summary_dir}/Metabrain_Augment_isoQTL_chr{chr}.txt',sep='\t')
isoQTL_lead=pd.read_csv(f'{qtl_lead_dir}/Metabrain_Augment_isoQTL_chr{chr}.txt',sep='\t')
isoQTL_sig_all=isoQTL_all[isoQTL_all['Probe'].isin(set(isoQTL_lead[isoQTL_lead['p']<5e-8]['Probe']))]

## read epiQTL
######## epiQTL ########
EpiQTL_list=[]
for qtl_type in epiQTltypes:
    QTL_sig=pd.read_csv(f'{epiQTL_sigProbe_summary_data}/{qtl_type}_sigProbe_summary_chr{chr}.txt',sep='\t')
    QTL_sig_all=QTL_sig[QTL_sig['SNP'].str.startswith('rs')]
    QTL_sig_all['qtl_type']=qtl_type
    EpiQTL_list.append(QTL_sig_all)
EpiQTL_sig_summary_df=pd.concat(EpiQTL_list,ignore_index=True)
EpiQTL_sig_summary_df['Probe']=EpiQTL_sig_summary_df['qtl_type']+'_'+EpiQTL_sig_summary_df['Probe']

df1=isoQTL_sig_all
df2=EpiQTL_sig_summary_df
common_snps = set(df1['SNP']) & set(df2['SNP'])
df1_overlap = df1[df1['SNP'].isin(common_snps)]
df2_overlap = df2[df2['SNP'].isin(common_snps)]
pair_counts = defaultdict(int)
df1_snp = df1_overlap.groupby("SNP")["Probe"].apply(list)
df2_snp = df2_overlap.groupby("SNP")["Probe"].apply(list)
common_snps = df1_snp.index.intersection(df2_snp.index)
for snp in common_snps:
    probes1 = df1_snp[snp]
    probes2 = df2_snp[snp]
    for p1, p2 in product(probes1, probes2):
        pair_counts[(p1, p2)] += 1
pair_df = pd.DataFrame([(p1, p2, n) for (p1, p2), n in pair_counts.items()],columns=["probe1","probe2","n_snp"])
pair_df = pair_df[pair_df["n_snp"] >= 100]
pair_df.columns=['isoform_probe','epi_probe','n_snp']
isoform_eff=df1_overlap[df1_overlap['Probe'].isin(set(pair_df['isoform_probe']))]
epi_dff=df2_overlap[df2_overlap['Probe'].isin(set(pair_df['epi_probe']))]
isoform_eff.to_csv(f'/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/isoQTL_epiQTL_input/overlap_summary_data/isoQTL_sig_Probe_summary_chr{chr}.txt',sep='\t',index=False)
epi_dff.to_csv(f'/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/isoQTL_epiQTL_input/overlap_summary_data/epiQTL_sig_Probe_summary_chr{chr}.txt',sep='\t',index=False)
pair_df.to_csv(f'/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/isoQTL_epiQTL_input/isoQTL_epiQTL_pair_Probe_chr{chr}',sep='\t',index='False')