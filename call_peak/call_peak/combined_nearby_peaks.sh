tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum bonemarrow ileum heart thymus stomach skin aorta tongue bladder CB jejunum uterus ovary BAT iWAT mammarygland)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
antibody=H3K27me3
for tissue in ${tissues[@]}
do
    bedtools merge -i ${data_path}/${tissue}/${antibody}/bed/${antibody}_young_old_merge-W1000-G3000-E100.bed -d 10000 | \
    awk 'BEGIN{OFS="\t"} {print $0, "peak" NR}' >  ${data_path}/${tissue}/${antibody}/bed/${antibody}_young_old_merge-W1000-G3000-E100_compress.bed
done