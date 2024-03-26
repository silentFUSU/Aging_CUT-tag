data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/

bam=($(ls ${data_path}*_ileum/bam/*.bam))  
featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_1kb_bins.saf -o ${data_path}/20240320_HJC_all_ileum/merge-1kb_bins.counts ${bam[@]} -F SAF -T 16