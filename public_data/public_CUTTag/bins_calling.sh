data_path=/mnt/transposon2/zhangyanxiaoLab/suzhuojie/project/Aging_CUT_Tag/public_data/werner_syndrom_aging/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
ref=mm10
files=$(ls ${data_path}bam/*.bam)
# featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/${ref}_10kb_bins.saf -o ${data_path}H3K9me3_10kb_bins.counts ${files} -F SAF -T 8 
featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/${ref}/${ref}_10kb_bins.saf -o ${data_path}H3K9me3_10kb_bins.counts ${files} -F SAF -T 8 
