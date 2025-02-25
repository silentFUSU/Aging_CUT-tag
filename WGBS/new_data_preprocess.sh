raw_data=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data_transposon2/20241211_combined_WGBS/
data_path=/mnt/transposon2/zhangyanxiaoLab/suzhuojie/project/Aging_CUT_Tag/samples/WGBS/
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/WGBS_search_table.csv
samples=$(ls ${raw_data}bed/*_CpG.bdg | sed 's|.*/||; s/_CpG\.bdg//') 
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
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/WGBS/WGBS_mkdir4samples.sh ${tissue}
    samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
    for sample_for_tissue in ${samples_for_tissue[@]}
    do
        if ls ${raw_data}bed/${sample_for_tissue}* 1> /dev/null 2>&1; then
            ln -s ${raw_data}bed/${sample_for_tissue}* ${data_path}${tissue}/bdg/
        fi

        if ls ${raw_data}bw/${sample_for_tissue}* 1> /dev/null 2>&1; then
            ln -s ${raw_data}bw/${sample_for_tissue}* ${data_path}${tissue}/bw/
        fi
    done
done