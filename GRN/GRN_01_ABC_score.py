### ABC score caculation

import os
import pandas as pd
import numpy as np
from multiprocessing import Pool

# extract tissue names
def get_tissues_from_matrix_dir(tissue_matrix_dir):
    tissue_files = [
        f for f in os.listdir(tissue_matrix_dir) 
        if f.endswith("_HiC.txt")
    ]
    tissues = [os.path.splitext(f)[0].replace("_HiC", "") for f in tissue_files]
    return tissues
    
# mapping peak location with bin location
def load_peak_bin_map(peak2bin_path):
    df = pd.read_csv(peak2bin_path, sep='\t', header=None, names=["bin", "peak"])
    return df.set_index("peak")["bin"].to_dict()

# information used in ABC score calculation
def load_abc_info(abc_info_path, peak2bin_map):
    df = pd.read_csv(abc_info_path, sep='\t', header=None, names=["gene", "tss", "peak", "dis"])
    df["dis"] = df["dis"].abs() 
    df["gene_bin"] = df["tss"].map(peak2bin_map)
    df["peak_bin"] = df["peak"].map(peak2bin_map)
    return df.dropna(subset=["gene_bin", "peak_bin"])


# Process ABC score calculation
def process_one_tissue(tissue_name, matrix_dir, activity_df, abc_info, peak2bin_map, out_path):
    hic_file = os.path.join(matrix_dir, f"{tissue_name}_HiC.txt")
    hic_data = pd.read_csv(hic_file, sep='\t', header=None, names=["bin1", "bin2", "count"], dtype={
        "bin1": np.int32, "bin2": np.int32, "count": np.float32
    })
    hic_data["bin_pair"] = hic_data["bin1"].astype(str) + "_" + hic_data["bin2"].astype(str)
    eg_df = abc_info.copy()
    eg_df["bin_pair"] = eg_df["gene_bin"].astype(str) + "_" + eg_df["peak_bin"].astype(str)
    eg_df["activity"] = eg_df["peak"].map(activity_df[tissue_name])
    pair2count = pd.Series(hic_data["count"].values, index=hic_data["bin_pair"].values)
    eg_df["contact"] = eg_df["bin_pair"].map(pair2count).fillna(0.0)
    eg_df["score"] = eg_df["activity"] * eg_df["contact"]
    tss_df = eg_df[["gene", "tss", "gene_bin"]].drop_duplicates()
    tss_df["activity"] = tss_df["tss"].map(activity_df[tissue_name])
    tss_df["contact"] = tss_df["gene_bin"].map(lambda b: pair2count.get(f"{b}_{b}", 0.0))
    tss_df["score"] = tss_df["activity"] * tss_df["contact"]
    tss_df["peak"] = "tss_dummy"
    eg_df = pd.concat([eg_df, tss_df], ignore_index=True)
    eg_df["denom"] = eg_df.groupby("gene")["score"].transform("sum").replace(0, np.nan)    
    eg_df["abc_score"] = eg_df["score"] / eg_df["denom"]
    eg_df = eg_df[eg_df['peak'] != "tss_dummy"]
    abc_result = eg_df[["gene", "peak", "abc_score"]].copy()
    abc_result.columns = ["gene", "peak", tissue_name]
    abc_result.to_csv(out_path, sep='\t', index=False)

# parallelism
def run_all(tissue_matrix_dir, activity_score_path, peak2bin_path, abc_info_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    activity_df = pd.read_csv(activity_score_path, sep='\t', index_col=0)
    peak2bin_map = load_peak_bin_map(peak2bin_path)
    abc_info = load_abc_info(abc_info_path, peak2bin_map)
    tissues = get_tissues_from_matrix_dir(tissue_matrix_dir)
    with Pool(processes=4) as pool:
        results = [
            pool.apply_async(process_one_tissue, args=(
                tissue,
                tissue_matrix_dir,
                activity_df,
                abc_info,
                peak2bin_map,
                os.path.join(out_dir, f"{tissue}_abc_score.txt")
            ))
            for tissue in tissues
        ]
        for r in results:
            r.get()

if __name__ == "__main__":
    run_all(
        tissue_matrix_dir="/path/to/HiC/file",
        #This folder is for HiC files named as "{tissue_name}_HiC.txt"
        activity_score_path="/path/to/abc_activity_score.txt",
        #This file is matrix of activity score of enhancers (enhancer x tissue)
        peak2bin_path="/path/to/peak2bin.txt",
        #This file is table for mapping enhancer with bin location
        abc_info_path="/path/to/abc_info.txt",
        #This file is table for information of each enhancer-gene pair, including gene name, TSS location, enhancer location, enhancer-gene distance
        out_dir="/path/to/ABC/score/result"
        #This folder is for ABC score results output
    )