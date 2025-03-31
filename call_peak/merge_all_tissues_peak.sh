tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum bonemarrow ileum heart thymus stomach skin aorta tongue bladder CB jejunum uterus ovary BAT iWAT mammarygland)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
antibody=H3K27me3
window_size=5000
gap_size=10000
e_value=100
files=()
for tissue in ${tissues[@]}
do
    bed=$(ls ${data_path}${tissue}/${antibody}/bed/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E100.bed)
    files+=("$bed")
done
cat  ${files[@]} | \
    sort -k1,1 -k2,2n | \
    bedtools merge  > ${data_path}all/${antibody}/bed/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}_bedtools.bed

awk '$3 - $2 > 100000' ${data_path}all/${antibody}/bed/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}_bedtools.bed > ${data_path}all/${antibody}/bed/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}_bedtools_filtered_100kb.bed  