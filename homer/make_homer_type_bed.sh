data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/ChromHMM/15_until_skin/state_transfer/E9_to_E1/bed/
beds=($(ls ${data_path}*.bed))
for bed in ${beds[@]}
do
    bed_name="${bed}"
    bed_new_name="${bed_name/.bed/_homer.bed}"
    awk '{print $1"\t"$2"\t"$3"\t""bin"NR"\t.\t."}' ${bed} > ${bed_new_name}
done