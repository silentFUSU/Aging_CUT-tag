tissues=(brain CB liver lung kidney colon)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
cleaned_file=$(mktemp)
cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file" 
for tissue in ${tissues[@]}
do
    samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
    for sample in ${samples_for_tissue[@]}
    do
        age=$(awk -F',' -v s="$sample" 'NR > 1 && ($3 == s) {print $5}' "$cleaned_file")  
        ln -s ${data_path}${tissue}/juicer/${sample}.allValidPairs.hic ${data_path}${tissue}/juicer/${sample}_${age}.allValidPairs.hic 
    done
done 

