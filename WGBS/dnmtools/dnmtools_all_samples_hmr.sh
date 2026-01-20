data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/
ref=/storage/zhangyanxiaoLab/xiongxiong/index/bismark/mm10/mm10.fa
tissue=$1

echo $tissue
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/WGBS_search_table.csv
cleaned_file=$(mktemp)  
cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
samples=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$cleaned_file")
IFS=$'\n' read -rd '' -a sample_array <<<"$samples"  
sample_counts=()
for sample in ${sample_array[@]}
do
    file=$(ls ${data_path}${tissue}/dnmtools_count/${sample}.counts)
    sample_counts+=("$file")
done
echo ${sample_counts[@]}
dnmtools hmr-rep ${sample_counts[@]} -o ${data_path}${tissue}/hmr/all_samples_hmr.bed
