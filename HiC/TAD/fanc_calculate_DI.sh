data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=lung
resolution=10kb
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
ref=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.fa
mkdir -p ${data_path}${tissue}/TAD/
mkdir -p ${data_path}${tissue}/TAD/fanc_DI
for sample in ${samples_for_tissue[@]}
do
    juicer_hic=${data_path}${tissue}/juicer/${sample}.allValidPairs.hic
    fanc directionality -o bed -w 2000000 -tmp \
        ${juicer_hic}@${resolution}@SCALE ${data_path}${tissue}/TAD/fanc_DI/${sample}_${resolution}.directionality &
done
wait 

for sample in ${samples_for_tissue[@]}
do
    bed=${data_path}${tissue}/TAD/fanc_DI/${sample}_${resolution}.directionality_2mb.bed
    awk -F'\t' '!(tolower($5) == "nan") {print $1 "\t" $2 "\t" $3 "\t" $5}' $bed > ${data_path}${tissue}/TAD/fanc_DI/${sample}_${resolution}.directionality_2mb.bedgraph
done
echo all done