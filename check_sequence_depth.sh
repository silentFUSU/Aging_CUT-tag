base_dir=$1
output_file=$2
samples=$(find "$base_dir" -type f -name "*R1.raw.fastq.gz")  
> "$output_file"  
for sample in ${samples[@]}
do
    line_count=$(zcat "$sample" | wc -l)  
    echo "$sample $line_count" >> "$output_file"  
done