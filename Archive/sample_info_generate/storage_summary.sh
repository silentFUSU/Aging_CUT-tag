data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/sample_info/ATAC_data_path.csv
cleaned_file=$(mktemp)
cat "$data_path" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
output_file=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/sample_info/ATAC_data_path_storage.csv

header=$(head -n 1 "$cleaned_file" )
echo -e "$header,V2_storage,V3_storage" > "$output_file"


tail -n +2 "$cleaned_file" | while IFS=',' read -r line; do
    # V2=$(echo "$line" | cut -d',' -f10  | sed -n 's/.*"\(.*\)".*/\1/p')
    V2=$(echo "$line" | cut -d',' -f11| sed -n 's/.*"\(.*\)".*/\1/p')
    # echo $V2
    # echo "$V2" | od -c
    # V3=$(echo "$line" | cut -d',' -f11 | sed -n 's/.*"\(.*\)".*/\1/p')
    V3=$(echo "$line" | cut -d',' -f12| sed -n 's/.*"\(.*\)".*/\1/p')

    convert_to_gb() {
        local size=$1
        if [[ $size == *M ]]; then
            # Remove 'M' and divide by 1000 to convert to gigabytes
            echo "scale=3; ${size%M}/1000" | bc
        elif [[ $size == *G ]]; then
            # Remove 'G' and keep the number as is
            echo "${size%G}"
        else
            echo "N/A" # Or handle other cases as needed
        fi
    }

    # Get space usage for V2 and V3
    V2_size=$(du -sh "$V2" 2>/dev/null | cut -f1)
    V3_size=$(du -sh "$V3" 2>/dev/null | cut -f1)

    V2_storage=$(convert_to_gb "$V2_size")
    V3_storage=$(convert_to_gb "$V3_size")

    echo "$line,$V2_storage,$V3_storage" >> "$output_file"
done