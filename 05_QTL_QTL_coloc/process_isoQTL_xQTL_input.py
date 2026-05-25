import pandas as pd
import os
from multiprocessing import Pool
import sys
chr = sys.argv[1]
out_dir='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/isoQTL_xQTL_input'
qtl_summary_dir = '/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_all'
qtl_lead_dir = '/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/Meta_analysis/Meta_QTL/query_per_probe_lead'
epiQTL_sigProbe_summary_data='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_unique_isoQTL_genetic_control/epi_sigQTL_data'
isoform_level_QTLtypes=['irQTL']
Gene_level_QTL_types=['RI_ratio','MX_ratio','SE_ratio','A3_ratio','A5_ratio','AF_ratio','AL_ratio','eQTL','sQTL']
isoQTL_all=pd.read_csv(f'{qtl_summary_dir}/Metabrain_Augment_isoQTL_chr{chr}.txt',sep='\t')
isoQTL_lead=pd.read_csv(f'{qtl_lead_dir}/Metabrain_Augment_isoQTL_chr{chr}.txt',sep='\t')
isoQTL_sig_all=isoQTL_all[isoQTL_all['Probe'].isin(set(isoQTL_lead[isoQTL_lead['p']<5e-8]['Probe']))]
irQTL_all=pd.read_csv(f'{qtl_summary_dir}/Metabrain_Augment_irQTL_chr{chr}.txt',sep='\t')
irQTL_lead=pd.read_csv(f'{qtl_lead_dir}/Metabrain_Augment_irQTL_chr{chr}.txt',sep='\t')
irQTL_sig_all=irQTL_all[irQTL_all['Probe'].isin(set(irQTL_lead[irQTL_lead['p']<5e-8]['Probe']))]
irQTL_sig_all['Probe']='irQTL'+'_'+irQTL_sig_all['Probe']
Gene_level_sig_qtl_summary_list=[]
for qtl_type in Gene_level_QTL_types:
    QTL_all=pd.read_csv(f'{qtl_summary_dir}/Metabrain_Augment_{qtl_type}_chr{chr}.txt',sep='\t')
    QTL_lead=pd.read_csv(f'{qtl_lead_dir}/Metabrain_Augment_{qtl_type}_chr{chr}.txt',sep='\t')
    if qtl_type=='sQTL':
        QTL_lead['Probe'] = QTL_lead['Probe'].str.replace('NOVELGENE', 'novelGene', regex=False)
        QTL_all['Probe'] = QTL_all['Probe'].str.replace('NOVELGENE', 'novelGene', regex=False)
        QTL_lead['Gene']= QTL_lead['Probe']
        QTL_all['Gene']= QTL_all['Probe']
    QTL_sig_all=QTL_all[QTL_all['Probe'].isin(set(QTL_lead[QTL_lead['p']<5e-8]['Probe']))]
    QTL_sig_all['qtl_type']=qtl_type
    Gene_level_sig_qtl_summary_list.append(QTL_sig_all)
Gene_level_sig_summary_df=pd.concat(Gene_level_sig_qtl_summary_list,ignore_index=True)
Gene_level_sig_summary_df['Probe']=Gene_level_sig_summary_df['qtl_type']+'_'+Gene_level_sig_summary_df['Probe']
Gene_level_sig_summary_df=Gene_level_sig_summary_df.drop(columns=['qtl_type'])

summary_list=[isoQTL_sig_all,irQTL_sig_all,Gene_level_sig_summary_df]
summary_sig_summary_df=pd.concat(summary_list,ignore_index=True)
summary_sig_summary_df=summary_sig_summary_df[summary_sig_summary_df['Gene'].isin(set(isoQTL_sig_all['Gene']))]
groups = summary_sig_summary_df.groupby('Gene')
def write_group(args):
    gene, df = args
    df.to_csv(f'{out_dir}/{gene}_sig_isoQTL_xQTL_summary.txt',
              sep='\t', index=False)
with Pool(32) as p:
    p.map(write_group, groups)