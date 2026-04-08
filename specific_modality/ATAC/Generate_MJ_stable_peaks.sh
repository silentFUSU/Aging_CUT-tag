Diff_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/ATAC_peak_from_LMJ/
all_peak_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/ATAC_peak_from_LMJ/tissue_set/
stable_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/ATAC_peak_from_LMJ/stable/
tissues=(skin aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas spleen stomach testis thymus tongue uterus mammarygland iWAT)
for tissue in ${tissues[@]}
do
    sort -k1,1 -k2,2n ${Diff_path}/up/${tissue}_Up.bed  > ${Diff_path}up/${tissue}_Up_sorted.bed
    sort -k1,1 -k2,2n ${Diff_path}/down/${tissue}_Down.bed > ${Diff_path}/down/${tissue}_Down_sorted.bed
    bedtools subtract -a ${all_peak_path}${tissue}.bed -b ${Diff_path}up/${tissue}_Up_sorted.bed  > ${stable_path}tmp.${tissue}_stable.bed
    bedtools subtract -a ${stable_path}tmp.${tissue}_stable.bed -b ${Diff_path}/down/${tissue}_Down_sorted.bed > ${stable_path}${tissue}_stable.bed
    rm ${stable_path}tmp.${tissue}_stable.bed
done