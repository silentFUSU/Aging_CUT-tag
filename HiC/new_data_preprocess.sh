raw_data=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data_transposon2/20250120_WJH_Thymus_HiC/
data_path=/mnt/transposon2/zhangyanxiaoLab/suzhuojie/project/Aging_CUT_Tag/samples/HiC/
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples=$(find ${raw_data}result/hic_results/data/ -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | grep -v tmp)  
cleaned_file=$(mktemp)  
cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
tissue_array=() 
while IFS= read -r prefix; do  
    tissue=$(awk -v prefix="$prefix" -F, '$3 == prefix {print $1}' "$cleaned_file")  
    if [ -n "$tissue" ]; then  
        tissue_array+=("$tissue")  
    fi  
done <<< "$samples" 
unique_tissue_array=($(printf "%s\n" "${tissue_array[@]}" | sort -u))  

for tissue in ${unique_tissue_array[@]}
do
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/HiC/HiC_mkdir4samples.sh ${tissue}
    samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
    for sample_for_tissue in ${samples_for_tissue[@]}
    do
        if ls ${raw_data}result/hic_results/data/${sample_for_tissue}/${sample_for_tissue}*.allValidPairs 1> /dev/null 2>&1; then  
            ln -s ${raw_data}result/hic_results/data/${sample_for_tissue}/${sample_for_tissue}*.allValidPairs ${data_path}${tissue}/ValidPairs/
        fi
        if ls ${raw_data}result/hic_results/data/${sample_for_tissue}*.allValidPairs.hic 1> /dev/null 2>&1; then  
            ln -s ${raw_data}result/hic_results/data/${sample_for_tissue}*.allValidPairs.hic ${data_path}${tissue}/juicer/
        fi
    done
done