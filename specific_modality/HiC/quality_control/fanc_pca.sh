tissue=lung
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/HiC/PCA/
resolution=100kb
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
cleaned_file=$(mktemp) 
cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
readarray -t young < <(awk -v tissue="$tissue" -F, '$1 == tissue && $5 == "3M" {print $3}' "$cleaned_file")  
readarray -t old < <(awk -v tissue="$tissue" -F, '$1 == tissue && $5 == "24M" {print $3}' "$cleaned_file")  

name_params=("${young[@]}" "${old[@]}")  
hic_files=()  
for sample in "${young[@]}"; do  
    hic_files+=("${data_path}${tissue}/juicer/${sample}.allValidPairs.hic@${resolution}")  
done  
for sample in "${old[@]}"; do  
    hic_files+=("${data_path}${tissue}/juicer/${sample}.allValidPairs.hic@${resolution}")  
done  

fanc pca -n  "${name_params[@]}" \
    -Z -s 1000000 -f -p ${result_path}${tissue}_FanC_${resolution}.pca.png \
    "${hic_files[@]}" ${result_path}${tissue}_${resolution}_FanC.pca

echo "fanc pca -n  "${name_params[@]}" -Z -s 1000000 -f -p ${result_path}${tissue}_FanC_${resolution}.pca.png "${hic_files[@]}" ${result_path}${tissue}_${resolution}_FanC.pca"
