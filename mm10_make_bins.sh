data_path="/storage/zhangyanxiaoLab/suzhuojie/ref_data/"
# bin_size=10000
# bin_size_label="10kb"
bin_size=1000
bin_size_label="1kb"

bedtools makewindows -g ${data_path}mm10_chromsize.bed -w ${bin_size} > ${data_path}tmp.mm10_${bin_size_label}_bins.bed
awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' ${data_path}tmp.mm10_${bin_size_label}_bins.bed |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}mm10_${bin_size_label}_bins.bed
bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}mm10_${bin_size_label}_bins.bed  ${data_path}mm10_${bin_size_label}_bins.saf
rm ${data_path}tmp.mm10_${bin_size_label}_bins.bed