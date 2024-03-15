data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
tissue=muscle
ref=mm10
antibody=H3K36me3
files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
