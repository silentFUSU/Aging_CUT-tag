black_list=/storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10-blacklist.v2.bed
bin_bed=/storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_1kb_bins.bed
bedtools intersect -a ${bin_bed} -b ${black_list} -wa > /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_1kb_bins_in_blacklist.bed
bin_bed=/storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_10kb_bins.bed
bedtools intersect -a ${bin_bed} -b ${black_list} -wa > /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_10kb_bins_in_blacklist.bed