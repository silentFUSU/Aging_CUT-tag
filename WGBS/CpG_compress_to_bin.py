import os  
import pandas as pd  
import numpy as np  
import multiprocessing as mp  
os.chdir("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")  

def process_chunk(chunk, ref):  
    ref_sorted = ref.sort_values(by=['V1', 'V2']).reset_index(drop=True)  
    chunk_sorted = chunk.sort_values(by=['V1', 'V2']).reset_index(drop=True)  
    merged = pd.merge_asof(chunk_sorted, ref_sorted, by='V1', left_on='V2', right_on='V2', direction='backward')  
    merged = merged[(merged['V2_y'] <= merged['V2_x']) & (merged['V3'] >= merged['V2_x'])]  
    chunk['bin_id'] = merged['V4']  
    return chunk  

def compress_to_bin(tissue, bin_size):  
    search_table = pd.read_csv("data/samples/all/WGBS_search_table.csv")  
    samples = search_table[search_table['tissue'] == tissue]['sample_name']  
    bin_methylation_summary = pd.DataFrame(columns=['bin_id', 'total_V4', 'total_V5', 'percent', 'sample'])  
    os.makedirs(f"data/samples/WGBS/{tissue}/compress2bin/", exist_ok=True)  
    ref = pd.read_csv(f"~/ref_data/mm10_{bin_size}_bins.bed", sep="\t", header=None, names=['V1', 'V2', 'V3', 'V4'])  
    ref = ref.sort_values(by=['V1', 'V2']).reset_index(drop=True)  
    for sample in samples:  
        df = pd.read_csv(f"data/samples/WGBS/{tissue}/bdg/{sample}_CpG.bdg", sep="\t", header=None, names=['V1', 'V2', 'V3', 'V4', 'V5', 'V6']) 
        # df['depth'] = df['V5']  
        chromosomes = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]
        df = df[df['V1'].isin(chromosomes)]  
        df = df.sort_values(by=['V1', 'V2']).reset_index(drop=True)  
        num_cores = 10  
        chunks = np.array_split(df, num_cores)  
        with mp.Pool(num_cores) as pool:  
            results = pool.starmap(process_chunk, [(chunk, ref) for chunk in chunks])  
        final_result = pd.concat(results, ignore_index=True)  
        final_result.to_csv(f"data/samples/WGBS/{tissue}/compress2bin/{sample}_CpG_compress_to_{bin_size}_bin_level_all_depth.csv", index=False)  
        bin_methylation = final_result.groupby('bin_id').agg({'V4': 'sum', 'V5': 'sum'}).reset_index()  
        bin_methylation['percent'] = bin_methylation['V4'] / bin_methylation['V5'] * 100  
        bin_methylation['sample'] = sample  
        bin_methylation_summary = pd.concat([bin_methylation_summary, bin_methylation], ignore_index=True)  
    
    bin_methylation_summary.to_csv(f"data/samples/WGBS/{tissue}/{bin_size}_bin_level_CpG_all_depth.csv", index=False)  
