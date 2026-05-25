import pandas as pd
import os
from concurrent.futures import ThreadPoolExecutor 
import re
import shutil
import argparse
### function for get sig result
def get_sig_result(QTL_dir,p_nominal=5e-8):
    work_dir='03_metabrain_cohort_data'
    basename=os.path.basename(QTL_dir)
    result_list=[]
    for i in range(1,23):
        res=pd.read_parquet(f'{QTL_dir}/{basename}.cis_qtl_pairs.{i}.parquet')
        res_sig=res[res['pval_nominal']<p_nominal]
        result_list.append(res_sig)
    result=pd.concat(result_list)
    return result

def run_get_sig_result(args):
    path, p_nominal = args
    return get_sig_result(path, p_nominal=p_nominal)

def run_sig_qtl(cohort,annotation,phenotype,p_nominal=5e-8, max_workers=10):
    work_dir=f'03_metabrain_cohort_data/{cohort}/Tensorqtl_out'
    PCA_work_dir = f'{work_dir}/{cohort}_{annotation}_{phenotype}_tier'
    parent_dir = os.path.dirname(PCA_work_dir)
    prefix_name = os.path.basename(PCA_work_dir)
    pc_pattern = re.compile(rf'{re.escape(prefix_name)}_PC(\d+)')
    PC_list = []
    for name in os.listdir(parent_dir):
        m = pc_pattern.match(name)
        if m:
            PC_list.append(int(m.group(1)))
    PC_list = sorted(PC_list)
    if not PC_list:
        raise ValueError(f'No PC directories/files found for prefix: {PCA_work_dir}')
    
    tasks = [(f'{PCA_work_dir}_PC{pc}', p_nominal) for pc in PC_list]
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(run_get_sig_result, tasks))
    # return dataframe dictionary
    return {f'PC{pc}': res for pc, res in zip(PC_list, results)}

def select_best_PC_for_maximizing_power(res_dict):
    Max_QTL_power=0  ## use sig SNP to measure the power
    Max_QTL_PC=0
    for PC, QTL_res in res_dict.items():
        QTL_power=len(QTL_res)
        if QTL_power>Max_QTL_power:
            Max_QTL_power=QTL_power
            Max_QTL_PC=PC
    return Max_QTL_PC, res_dict[Max_QTL_PC]

def main(cohort,annotation,phenotype,thread,p_cutoff):
    work_dir='03_metabrain_cohort_data'
    tier_QTL_res=run_sig_qtl(cohort,annotation,phenotype,p_nominal=p_cutoff,max_workers=thread)
    best_pc,best_QTL_res=select_best_PC_for_maximizing_power(tier_QTL_res)
    os.makedirs(f"{work_dir}/{cohort}/Tensorqtl_best", exist_ok=True)
    os.makedirs(f"{work_dir}/{cohort}/Tensorqtl_best/{cohort}_{annotation}_{phenotype}", exist_ok=True)
    out_put_dir=f'{work_dir}/{cohort}/Tensorqtl_best/{cohort}_{annotation}_{phenotype}'
    ## save sig result
    best_QTL_res.to_csv(f'{out_put_dir}/{cohort}_{annotation}_{phenotype}_{best_pc}.sig.txt',sep='\t',index=False)
    ## mv snp-phenotype paired result
    shutil.move(f"{work_dir}/{cohort}/Tensorqtl_out/{cohort}_{annotation}_{phenotype}_tier_{best_pc}", f"{out_put_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert sample names to individual names and make bed file")
    parser.add_argument("--cohort",required=True,help="cohort name(CMC,ROSMAP et al)")
    parser.add_argument("--annotation",required=True,help="Gencode,Metabrain,Metabrain_Augment")
    parser.add_argument("--phenotype",required=True,help="such as gene count, transcript count,transcript ratio")
    parser.add_argument("--thread",default=10,help="threads used to parallerl")
    parser.add_argument("--Pcutoff",default=5e-8,help="nominal P value cutoff to define the signifcaiton association between variant and phenotype") 
    args = parser.parse_args()
    main(cohort=args.cohort,annotation=args.annotation,phenotype=args.phenotype,thread=args.thread,p_cutoff=args.Pcutoff)

### usage ###
# python script.py --cohort --annotation --phenotype --thread --Pcutoff