Diff_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/ATAC_peak_from_MJ/Diff_peaks/
all_peak_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/ATAC_peak_from_MJ/all_peaks/
stable_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/ATAC_peak_from_MJ/stable_peaks/
tissues=(skin aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas spleen stomach testis thymus tongue uterus mammarygland iWAT)
for tissue in ${tissues[@]}
do
    sort -k1,1 -k2,2n ${Diff_path}Peak_Up_${tissue}.bed  > ${Diff_path}Peak_Up_${tissue}_sorted.bed 
    sort -k1,1 -k2,2n ${Diff_path}Peak_Down_${tissue}.bed > ${Diff_path}Peak_Down_${tissue}_sorted.bed
    bedtools subtract -a ${all_peak_path}${tissue}*.bed -b ${Diff_path}Peak_Up_${tissue}_sorted.bed  > ${stable_path}tmp.${tissue}_stable.bed
    bedtools subtract -a ${stable_path}tmp.${tissue}_stable.bed -b ${Diff_path}Peak_Down_${tissue}_sorted.bed > ${stable_path}${tissue}_stable_peak.bed
    rm ${stable_path}tmp.${tissue}_stable.bed
done