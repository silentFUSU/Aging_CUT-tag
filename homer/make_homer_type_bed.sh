data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/ATAC_peak_from_MJ/
beds=($(ls ${data_path}*.bed))
for bed in ${beds[@]}
do
    bed_name="${bed}"
    bed_new_name="${bed_name/.bed/_homer.bed}"
    awk '{print $1"\t"$2"\t"$3"\t""peak"NR"\t.\t."}' ${bed} > ${bed_new_name}
done


data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/ATAC/motif_bg/bed/
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
for tissue in ${tissues[@]}
do 
    ln -s ${data_path}${tissue}/ATAC/bed/ATAC_peaks_diff_up.bed ${result_path}${tissue}_ATAC_peaks_diff_up.bed
    ln -s ${data_path}${tissue}/ATAC/bed/ATAC_peaks_diff_down.bed ${result_path}${tissue}_ATAC_peaks_diff_down.bed
    ln -s ${data_path}${tissue}/ATAC/bed/ATAC_peaks_diff_stable.bed ${result_path}${tissue}_ATAC_peaks_diff_stable.bed
done