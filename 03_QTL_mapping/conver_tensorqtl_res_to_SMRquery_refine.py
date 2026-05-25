### make transcript expression query format
from multiprocessing import process
import os
import glob
from selectors import EpollSelector
import pandas as pd
import numpy as np
import argparse
output_dir='cohort_query_format'
work_dir='03_metabrain_cohort_data'
transcript_info_dir='transcript_info_chr'
gene_info_dir='gene_info_chr'
cohorts=['Mayo','LIBD','HBCC','ROSMAP','CMC','BrainGVEx','BM10','BM22','BM36','BM44']
def main(chr,phenotype):
    print(f'process chr{chr}')
    bad_SNP=pd.read_csv('bad_rsids.txt',sep='\t',names=['SNP'])
    transcript_info=pd.read_csv(f'{transcript_info_dir}/transcript_chr{chr}.txt',sep='\t',usecols=[0,3,4,5,6,7])
    gene_info=pd.read_csv(f'{gene_info_dir}/gene_chr{chr}.txt',sep='\t')
    for cohort in cohorts:
        #cohort='ROSMAP'
        bim_file=f'{work_dir}/{cohort}/genotype_hg38/{cohort}_chr{chr}.bim'
        bim=pd.read_csv(bim_file,sep='\t',usecols=[0,1,3,4,5],names=['chr','variant_id','position','A1','A2'])
        for annotation in ['Gencode','Metabrain_Augment']:
            if phenotype in ['transcript_count','transcript_ratio']:
                QTL_res_dir=f'{work_dir}/{cohort}/Tensorqtl_best/{cohort}_{annotation}_{phenotype}'
                pc_dirs = glob.glob(os.path.join(QTL_res_dir, f"{cohort}_{annotation}_{phenotype}_tier_PC*"))
                best_pc_dir = pc_dirs[0]
                qtl_res_prefix = os.path.basename(best_pc_dir)
                qtl_res=pd.read_parquet(f'{best_pc_dir}/{qtl_res_prefix}.cis_qtl_pairs.{chr}.parquet')
                qtl_res_bim=pd.merge(qtl_res,bim,how='inner',on='variant_id')
                qtl_res_bim_transcript=pd.merge(qtl_res_bim,transcript_info,how='inner',left_on='phenotype_id',right_on='transcript_id')
                qtl_res_bim_transcript_use=qtl_res_bim_transcript[['variant_id','chr','position','A1','A2','af','phenotype_id','chr','TSS','gene_id','strand','slope','slope_se','pval_nominal']]
                qtl_res_bim_transcript_use.columns=['SNP','Chr','BP','A1','A2','Freq','Probe','Probe_Chr','Probe_bp','Gene','Orientation','b','se','p']
                qtl_res_bim_transcript_use=qtl_res_bim_transcript_use[~qtl_res_bim_transcript_use['SNP'].isin(bad_SNP['SNP'])]
                if phenotype=='transcript_count':
                    QTL_type='isoQTL'
                if phenotype=='transcript_ratio':
                    QTL_type='irQTL'
                qtl_res_bim_transcript_use.dropna().to_csv(f'{output_dir}/{annotation}/{cohort}_{QTL_type}_chr{chr}.txt',sep='\t',index=False,float_format='%.6g')
            
            if phenotype in ['gene_count','A3_ratio','A5_ratio','AL_ratio','AF_ratio','SE_ratio','RI_ratio','MX_ratio']:
                QTL_res_dir=f'{work_dir}/{cohort}/Tensorqtl_best/{cohort}_{annotation}_{phenotype}'
                pc_dirs = glob.glob(os.path.join(QTL_res_dir, f"{cohort}_{annotation}_{phenotype}_tier_PC*"))
                best_pc_dir = pc_dirs[0]
                qtl_res_prefix = os.path.basename(best_pc_dir)
                qtl_res=pd.read_parquet(f'{best_pc_dir}/{qtl_res_prefix}.cis_qtl_pairs.{chr}.parquet')
                qtl_res_bim=pd.merge(qtl_res,bim,how='inner',on='variant_id')
                qtl_res_bim['gene_id'] = qtl_res_bim['phenotype_id'].str.split(';').str[0]
                qtl_res_bim_gene=pd.merge(qtl_res_bim,gene_info,how='inner',on='gene_id')
                qtl_res_bim_gene_use=qtl_res_bim_gene[['variant_id','chr','position','A1','A2','af','phenotype_id','chr','gene_TSS','gene_id','strand','slope','slope_se','pval_nominal']]
                qtl_res_bim_gene_use.columns=['SNP','Chr','BP','A1','A2','Freq','Probe','Probe_Chr','Probe_bp','Gene','Orientation','b','se','p']
                ## remove bed snp with uncertainty position
                qtl_res_bim_gene_use=qtl_res_bim_gene_use[~qtl_res_bim_gene_use['SNP'].isin(bad_SNP['SNP'])] 
                ## remove association with empety value
                if phenotype=='gene_count':
                    QTL_type='eQTL'
                else:
                    QTL_type=phenotype
                qtl_res_bim_gene_use.dropna().to_csv(f'{output_dir}/{annotation}/{cohort}_{QTL_type}_chr{chr}.txt',sep='\t',index=False,float_format='%.6g')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert sample names to individual names and make bed file")
    parser.add_argument("--chr",required=True,help="1-22")
    parser.add_argument("--phenotype",required=True,help="transcript_count,transcript_ratio,gene_count,A5_ratio......") 
    args = parser.parse_args()
    main(chr=args.chr,phenotype=args.phenotype)