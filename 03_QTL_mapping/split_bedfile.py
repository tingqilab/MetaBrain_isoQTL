import os
import pandas as pd
from pathlib import Path
import argparse
import numpy as np
def split_bed_by_chr(folder_path):
    folder_path = Path(folder_path)
    bed_files = list(folder_path.glob("*.bed"))  # 获取所有 bed 文件
    for bed_file in bed_files:
        # 子文件夹名 = bed 文件名去掉 .bed
        subfolder = folder_path / bed_file.stem
        subfolder.mkdir(parents=True, exist_ok=True)
        # 读取 bed 文件
        bed_df = pd.read_csv(bed_file, sep='\t')
        # 确保 chrom 列存在，并重命名为 #chr
        if 'chrom' in bed_df.columns:
            bed_df.rename(columns={'chrom': '#chr'}, inplace=True)
        # 按染色体拆分
        for i in range(1, 23):
            chr_df = bed_df[bed_df['#chr'] == i]
            chr_file = subfolder / f"chr{i}.bed"
            chr_df.to_csv(chr_file, sep='\t', index=False)
    print(f"Processed {len(bed_files)} bed files in {folder_path}")

def main(Residual_dir):
    print(f"start process for cohort {Residual_dir}")
    split_bed_by_chr(Residual_dir)
    print(f"ending process for cohort {Residual_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="split phenotype Residuals bed file into 22 chrom")
    parser.add_argument("--Residual",required=True,help="Residual dir which include many bed file ")
    args = parser.parse_args()
    main(args.Residual)


