### split epigenomic isoQTL and isoQTL COLOC Data
import pandas as pd
import sys
chr = sys.argv[1]
work_dir='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/isoQTL_epiQTL_input'
paired_probe=pd.read_csv(f'{work_dir}/pair_Probe/isoQTL_epiQTL_pair_Probe_chr{chr}',sep='\t')
paired_probe=paired_probe.drop(columns=['Unnamed: 0'])
isoQTL_summary=pd.read_csv(f'{work_dir}/overlap_summary_data/isoQTL_sig_Probe_summary_chr{chr}.txt',sep='\t')
epiQTL_summary=pd.read_csv(f'{work_dir}/overlap_summary_data/epiQTL_sig_Probe_summary_chr{chr}.txt',sep='\t')
split_probe_pair_dir='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/isoQTL_epiQTL_input/split_pair_Probe'
split_summary_dir='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/05_trait_association/COLOC/isoQTL_epiQTL_input/split_overlap_summary_data'
chunk_size = 100000
paired_probe_chunks = [
    paired_probe.iloc[i:i+chunk_size]
    for i in range(0, len(paired_probe), chunk_size)
]
for index, pair in enumerate(paired_probe_chunks, start=1):
    subset_isoQTL = isoQTL_summary[isoQTL_summary['Probe'].isin(set(pair['isoform_probe']))]
    subset_epiQTL = epiQTL_summary[epiQTL_summary['Probe'].isin(set(pair['epi_probe']))]
    subset_isoQTL.to_csv(f'{split_summary_dir}/isoQTL_summary_chr{chr}_{index}.txt', sep='\t', index=False)
    subset_epiQTL.to_csv(f'{split_summary_dir}/epiQTL_summary_chr{chr}_{index}.txt', sep='\t', index=False)
    pair.to_csv(f'{split_probe_pair_dir}/Probe_pair_chr{chr}_{index}.txt', sep='\t')