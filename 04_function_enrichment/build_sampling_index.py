########## 1. index all QTL probe-MAF_bin-TSS_bin #################
import pandas as pd
import argparse
def main(QTL_type,chr):
    meta_qtl_dir='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/03_metabrain_cohort_data/Meta_analysis/Meta_QTL'
    index_out_dir='/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_function_enrichment/index_variant1'
    sig_probe=pd.read_csv(f'{meta_qtl_dir}/query_top/Metabrain_Augment_{QTL_type}_chr{chr}.cis.summary.txt',sep='\t',usecols=[0])
    sig_probe_set = set(sig_probe['Probe'])
    all_qtl=pd.read_csv(f'{meta_qtl_dir}/query_all/Metabrain_Augment_{QTL_type}_chr{chr}.txt',sep='\t')
    all_qtl=all_qtl[all_qtl['Probe'].isin(sig_probe_set)] ### only including the sig Probe
    all_qtl['TSS_dis'] = abs(all_qtl['Probe_bp'] - all_qtl['BP'])
    maf_bins = [0, 0.01, 0.05, 0.10, 0.20, 0.30, 0.50]
    tss_dis_bins=[0,1000,10000,100000,1000000]
    maf_labels = [f"{maf_bins[i]}_{maf_bins[i+1]}" for i in range(len(maf_bins)-1)]
    tss_labels = [f"{tss_dis_bins[i]}_{tss_dis_bins[i+1]}" for i in range(len(tss_dis_bins)-1)]
    all_qtl['TSS_bin']=pd.cut(all_qtl['TSS_dis'], bins=tss_dis_bins, labels=tss_labels, include_lowest=True)
    all_qtl['MAF_bin'] = pd.cut(all_qtl['Freq'], bins=maf_bins, labels=maf_labels, include_lowest=True)
    all_qtl['index'] = (all_qtl['TSS_bin'].astype(str)+ '-' + all_qtl['MAF_bin'].astype(str)+ '-' + all_qtl['Probe'].astype(str))
    all_qtl.to_csv(f'{index_out_dir}/{QTL_type}_chr{chr}_index.txt',sep='\t',index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="bulid QTL index")
    parser.add_argument("--QTLtype",required=True,help="isoQTL,eQTL,irQTL,sQTL")
    parser.add_argument("--chr",required=True,help="1...22") 
    args = parser.parse_args()
    main(QTL_type=args.QTLtype,chr=args.chr)