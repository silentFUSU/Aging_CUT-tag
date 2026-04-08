# directory=/storage/zhangyanxiaoLab/fastq/2023/2023-12-08-Jiangbei-LLX/
# output_file=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/BIG_list/md5_results.txt
# > "$output_file"
# find "$directory" -type f -name "*.gz" | while read -r file; do
#     md5_value=$(md5sum "$file" | awk '{ print $1 }')
#     file_name=$(basename "$file")
#     echo -e "$file_name\t$md5_value" >> "$output_file"
# done

data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/sample_info/HiC_data_path.csv
cleaned_file=$(mktemp)
cat "$data_path" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
output_file=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/sample_info/HiC_data_path_md5.csv

header=$(head -n 1 "$cleaned_file" )
# echo -e "$header,V2_md5,V3_md5" > "$output_file"
echo -e "$header" > "$output_file"

tail -n +2 "$cleaned_file" | while IFS=',' read -r line; do
    V2=$(echo "$line" | cut -d',' -f10  | sed -n 's/.*"\(.*\)".*/\1/p')
    # echo "$V2" | od -c
    V3=$(echo "$line" | cut -d',' -f11 | sed -n 's/.*"\(.*\)".*/\1/p')
    # echo "$V3" | od -c
    V2_md5=""
    V3_md5=""
    if [[ -f "$V2" ]]; then
        V2_md5=$(md5sum "$V2" | awk '{print $1}')
        echo $V2_md5
    fi
    if [[ -f "$V3" ]]; then
        V3_md5=$(md5sum "$V3" | awk '{print $1}')
        echo $V3_md5
    fi
    echo "$line,$V2_md5,$V3_md5" >> "$output_file"
done