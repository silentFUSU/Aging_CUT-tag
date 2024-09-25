tissue=kidney
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/${tissue}/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/WGBS/${tissue}/
CUTTag_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/

bed=${CUTTag_path}${tissue}/H3K27me3_H3K9me3_intersect/10kb_all_significant_second_quadrant.bed
conditions=(increase decrease)
for condition in ${conditions[@]}
do
    WGBS_bed=${data_path}DSS_table/bed/${tissue}_DMR_${condition}.bed
    bedtools intersect -a ${WGBS_bed} \
            -b ${bed} -wa > ${data_path}DSS_table/bed/${tissue}_DMR_${condition}_overlap_heterochromatin_switch.bed 
done