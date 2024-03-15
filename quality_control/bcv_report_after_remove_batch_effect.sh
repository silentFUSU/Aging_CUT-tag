data_path=$1
tissue=$2
# /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
antibodys=(H3K27me3 H3K9me3 H3K36me3 H3K27ac H3K4me3 H3K4me1)
output_file=${data_path}${tissue}/bcv_report_after_remove_batch_effect.txt
> ${output_file}

for antibody in ${antibodys[@]}
do
    file_path=$(ls ${data_path}${tissue}/${antibody}/${antibody}_{macs_narrowpeak,merge-W1000-G3000-E100}_diff_after_remove_batch_effect.csv)
    tail -n +2 ${file_path} | awk -v sample="$antibody"  -F ',' '{print $12, sample}' >> $output_file
done