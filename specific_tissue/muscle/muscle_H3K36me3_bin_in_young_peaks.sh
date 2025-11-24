data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
tissue=muscle
ref=mm10
antibody=H3K36me3
bam=$(ls ${data_path}${tissue}/${antibody}/bam/{LLX286,LLX298}*.bam)
age=old
samtools merge -o ${data_path}${tissue}/${antibody}/tmp.${age}_merge.bam ${bam} -@ 16
samtools index ${data_path}${tissue}/${antibody}/tmp.${age}_merge.bam -@ 16

window_size=1000
gap_size=3000
e_value=100
sicer  -t ${data_path}${tissue}/${antibody}/tmp.${age}_merge.bam  -o ${data_path}${tissue}/${antibody}/peaks \
    -s ${ref} -w ${window_size} -rt 16 -f 300 -egf 0.8 -fdr 0.01 -g ${gap_size} -e ${e_value} -cpu 21

awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' ${data_path}${tissue}/${antibody}/peaks/tmp.${age}_merge-W${window_size}-G${gap_size}.scoreisland \
    |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/bed/${antibody}_${age}_merge-W${window_size}-G${gap_size}-E${e_value}.bed

bin_size=10kb
bedtools intersect -a ${data_path}${tissue}/${antibody}/bed/${antibody}_${age}_merge-W${window_size}-G${gap_size}-E${e_value}.bed \
    -b ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_up.bed -wa > ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_up_in_${age}_peaks.bed

python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/H3K27me3_H3K9me3_relationship/remove_duplicate.py ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_up_in_${age}_peaks.bed

sort -k 1,1 -k2,2n ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_up_in_${age}_peaks_unique.bed >  ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_up_in_${age}_peaks_unique_sort.bed

rm ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_up_in_${age}_peaks.bed
rm ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_up_in_${age}_peaks_unique.bed

bedtools intersect -a ${data_path}${tissue}/${antibody}/bed/${antibody}_${age}_merge-W${window_size}-G${gap_size}-E${e_value}.bed \
    -b ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_down.bed -wa > ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_down_in_${age}_peaks.bed

python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/H3K27me3_H3K9me3_relationship/remove_duplicate.py ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_down_in_${age}_peaks.bed

sort -k 1,1 -k2,2n ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_down_in_${age}_peaks_unique.bed >  ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_down_in_${age}_peaks_unique_sort.bed

rm ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_down_in_${age}_peaks.bed
rm ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_down_in_${age}_peaks_unique.bed

rm ${data_path}${tissue}/${antibody}/tmp*