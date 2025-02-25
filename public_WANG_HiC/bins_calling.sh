data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/cellular_aging_GSE133292/Chip_seq/
samples=(SRR13274568 SRR13274569 SRR13274574 SRR13274575) 
files=$(ls ${data_path}bam/{SRR13274568,SRR13274569,SRR13274574,SRR13274575}.nodup.bam)
featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/hg38/hg38_10kb_bins.saf -o ${data_path}H3K27me3_10kb_bins.counts ${files} -F SAF -T 8 