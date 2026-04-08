data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
# tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum bonemarrow ileum heart thymus)
tissues=(skin stomach)
# tissue=bonemarrow
ref=mm10
quadrants=(first second third fourth)
board_bin_size=10kb
peak_bin_size=1kb
antibodys=(H3K27ac H3K4me1 H3K4me3)
for antibody in ${antibodys[@]}
do
    for tissue in ${tissues[@]}
    do
        for quadrant in ${quadrants[@]}
        do
        bedtools intersect -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_${peak_bin_size}_bins.bed \
            -b ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/${board_bin_size}_all_significant_${quadrant}_quadrant.bed -wa \
            >  ${data_path}${tissue}/${antibody}/bed/tmp_H3K27me3_H3K9me3_all_significant_${peak_bin_size}_${quadrant}_quadrant.bed

        python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/H3K27me3_H3K9me3_relationship/remove_duplicate.py \
            ${data_path}${tissue}/${antibody}/bed/tmp_H3K27me3_H3K9me3_all_significant_${peak_bin_size}_${quadrant}_quadrant.bed

        sort -k 1,1 -k2,2n ${data_path}${tissue}/${antibody}/bed/tmp_H3K27me3_H3K9me3_all_significant_${peak_bin_size}_${quadrant}_quadrant_unique.bed \
            >  ${data_path}${tissue}/${antibody}/bed/H3K27me3_H3K9me3_all_significant_${peak_bin_size}_${quadrant}_quadrant_unique_sort.bed

        rm ${data_path}${tissue}/${antibody}/bed/tmp*
        done
    done
done