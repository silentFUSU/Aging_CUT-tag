ln -s /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/lung/H3K9me3/bam/HJC*.bam* /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data_transposon2/20250110_DYQ_CUTTag/bam/

data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data_transposon2/20250110_DYQ_CUTTag/
ref=mm10
files=$(ls ${data_path}bam/*.bam)
featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_10kb_bins.saf -o ${data_path}H3K9me3_10kb_bins.counts ${files} -F SAF -T 8 