#### generate the annot file of obs snp and null snp for SLDSC, it must match baseline order
import pandas as pd
import argparse
from pathlib import Path

# ========== 参数解析 ==========
parser = argparse.ArgumentParser(description="Generate annotation files for one chromosome.")
parser.add_argument("--chr", type=int, required=True, help="Chromosome number (1-22)")
parser.add_argument("--QTL_type", type=str, required=True, help="Type of QTL") ## Brain_meta_sQTL_Augment, Brain_meta_eQTL_Augment,Brain_meta_sQTL_Gencode, Brain_meta_eQTL_Gencode
args = parser.parse_args()
chr_num = args.chr
QTL_type = args.QTL_type

# ========== 路径设置 ==========
baseline_v1_prefix = '/qitinglab/sharedata/LDSC_resource/GRCh38/baseline_v1.2/baseline.'
null_annot_dir = '/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_trait_hertability_enrichment/SLDSC/null_annot'
obs_annot_dir = '/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_trait_hertability_enrichment/SLDSC/obs_annot'
sampling_res= f'/qitinglab/shenlongxin/metabrain/PCA_based_qtl_mapping/04_trait_hertability_enrichment/SLDSC/sampling_res/{QTL_type}_sampling.txt'

# ========== 读取数据 ==========
print(f"==> Start processing chr{chr_num}")
sampling_res=pd.read_csv(sampling_res,sep='\t')
# baseline
baseline_df = pd.read_csv(f'{baseline_v1_prefix}{chr_num}.annot.gz', sep='\t', usecols=[0,1,2,3])

# obs annotation
# baseline_df[f'{QTL_type}'] = baseline_df['SNP'].isin(sampling_res['lead_SNP']).astype(int)
# obs_out = f'{obs_annot_dir}/{QTL_type}_obs.{chr_num}.annot.gz'
# baseline_df.to_csv(obs_out, sep='\t', index=False, compression='gzip')
# print(f"Saved obs annotation for chr{chr_num}: {obs_out}")

# null annotation
print(f"Start generating null annotations for chr{chr_num} ...")
# 提前把每列转成 set，加快 isin 判断
for rep in range(1, 101):
    snp_set = set(sampling_res[f'iter_{rep}'])
    baseline_df[f'{QTL_type}'] = baseline_df['SNP'].isin(snp_set).astype(int)
    out_file = f'{null_annot_dir}/{QTL_type}_rep{rep}.{chr_num}.annot.gz'
    baseline_df.to_csv(out_file, sep='\t', index=False, compression='gzip')
    if rep % 50 == 0:
        print(f"  → saved {rep}/1000 null files")
print(f"==> Finished chr{chr_num}")
