import snapatac2 as snap
import re
from pathlib import Path
import pandas as pd
import os
import polars as pl
import argparse

def read_bed_file(bed_file):
    regions_dict = {}
    bed_file = Path(bed_file)
    bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chr", "start", "end"])
    bed_df["region"] = bed_df["chr"] + ":" + bed_df["start"].astype(str) + "-" + bed_df["end"].astype(str)
    regions_dict[bed_file.name] = bed_df["region"].tolist()
    return regions_dict

def save_enrichment_results(enrichment_results, output_folder):
    for group_name, df in enrichment_results.items():
        output_file = os.path.join(output_folder, f"enrichment_results_{group_name}.csv")
        df_pandas = df.to_pandas()
        df_pandas.to_csv(output_file, index=True)

genome_fasta="/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.fa"
method = "binomial"
# motifs = snap._snapatac2.read_motifs("/storage/zhangyanxiaoLab/suzhuojie/ref_data/TE_reference/motif_databases/CIS-BP_2.00/Mus_musculus.meme")
motifs = snap._snapatac2.read_motifs("/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/meme_cisbp.meme")
result_path = "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/ATAC_peak_from_LMJ/motif_snapatac2_cisbp/motif_bg_db_from_PFM/"

parser = argparse.ArgumentParser()
parser.add_argument('--condition', type=str, required=True, help='Specify the condition')
parser.add_argument('--tissue', type=str, required=True, help='Specify the tissue type')
args = parser.parse_args()
tissue = args.tissue
condition = args.condition
print(tissue)
print(condition)
regions = read_bed_file(f"/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/ATAC_peak_from_LMJ/{condition}/{tissue}_{condition.capitalize()}_sorted.bed")
background_bed = f"/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/ATAC_peak_from_LMJ/stable/{tissue}_stable.bed"
# background_bed = f"/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/ATAC_peak_from_LMJ/tissue_set/{tissue}.bed"
bed_df_back = pd.read_csv(background_bed , sep="\t", header=None, names=["chr", "start", "end"])
bed_df_back["region"] = bed_df_back["chr"] + ":" + bed_df_back["start"].astype(str) + "-" + bed_df_back["end"].astype(str)
background = bed_df_back["region"].tolist()
enrichment_results = snap.tl.motif_enrichment(motifs, regions, genome_fasta, background=background, method=method)
save_enrichment_results(enrichment_results, f"{result_path}{condition}/")