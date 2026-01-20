directory=/storage/zhangyanxiaoLab/fastq/2023/2023-12-08-Jiangbei-LLX/
output_file=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/BIG_list/md5_results.txt
> "$output_file"
find "$directory" -type f -name "*.gz" | while read -r file; do
    md5_value=$(md5sum "$file" | awk '{ print $1 }')
    file_name=$(basename "$file")
    echo -e "$file_name\t$md5_value" >> "$output_file"
done