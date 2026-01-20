data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/
ref=/storage/zhangyanxiaoLab/xiongxiong/index/bismark/mm10/mm10.fa
tissue=$1

echo $tissue
mkdir -p ${data_path}${tissue}/dnmtools_count/
mkdir -p ${data_path}${tissue}/hmr/
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/WGBS_search_table.csv
cleaned_file=$(mktemp)  
cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
samples=$(awk -F',' -v t="$tissue" -v y="3M" 'NR > 1 && ($1 == t) && ($5 == y) {print $3}' "$cleaned_file")
IFS=$'\n' read -rd '' -a sample_array <<<"$samples"  
echo ${sample_array[@]}
for sample in ${sample_array[@]}
do
    dnmtools counts -c ${ref} ${data_path}${tissue}/bam/${sample}*_dedup.bam -o ${data_path}${tissue}/dnmtools_count/${sample}.counts &
done
wait
young=()
for sample in ${sample_array[@]}
do
    file=$(ls ${data_path}${tissue}/dnmtools_count/${sample}.counts)
    young+=("$file")
done
echo ${young[@]}
dnmtools hmr-rep ${young} -o ${data_path}${tissue}/hmr/young_hmr.bed
