data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/
beds=($(ls ${data_path}*.bed))
for bed in ${beds[@]}
do
    bed_name="${bed}"
    bed_new_name="${bed_name/.bed/_homer.bed}"
    awk '{print $1"\t"$2"\t"$3"\t""peak"NR"\t.\t."}' ${bed} > ${bed_new_name}
done

data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/H3K9me3/recursion_peaks_diff_table/
mkdir ${data_path}motif_bg

kmeans=(1 2 3 4)
for kmean in ${kmeans[@]}
do
    input=${data_path}bed/kmeans${kmean}_uinon_recursion_peaks_homer.bed
    find "${data_path}" -type f -name '*_peaks.bed' \
        ! -path "${data_path}bed/kmeans${kmean}_uinon_recursion_peaks.bed" \
        -exec cat {} + |
        sort -k1,1 -k2,2n > ${data_path}bed/tmp_merged_sorted_bed_file.bed
    bed=${data_path}bed/tmp_merged_sorted_bed_file.bed
    bed_name="${data_path}bed/tmp_merged_sorted_bed_file.bed"
    bed_new_name="${bed_name/.bed/_homer.bed}"
    awk '{print $1"\t"$2"\t"$3"\t""peak"NR"\t.\t."}' ${bed} > ${bed_new_name}
    
    background=${data_path}bed/tmp_merged_sorted_bed_file_homer.bed
    findMotifsGenome.pl ${input} mm10 ${data_path}motif_bg/kmeans${kmean} -size 200 -mask -bg ${background} 
done

mkdir ${data_path}motif_without_bg
kmeans=(1 2 3 4)
for kmean in ${kmeans[@]}
do
    input=${data_path}bed/kmeans${kmean}_uinon_recursion_peaks_homer.bed
    findMotifsGenome.pl ${input} mm10  ${data_path}motif_without_bg/kmeans${kmean} -size 200 -mask 
done