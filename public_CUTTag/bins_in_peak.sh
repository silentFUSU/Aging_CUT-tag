data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/Hippocampus_aging/

window_size=1000
gap_size=3000
e_value=100
ref_data=~/ref_data/
bin_size=10kb
antibody=H3K9me3
bedtools intersect -a ${ref_data}mm10_${bin_size}_bins.bed \
    -b ${data_path}/bed/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed -wa > ${data_path}/bed/${antibody}_${bin_size}_in_young_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed


