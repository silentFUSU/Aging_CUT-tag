add_data_path=/storage/zhangyanxiaoLab/fastq/2024/2024-10-26-Jiangbei-DYQ/
previous_data_path_list=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/WGBS_data_path.csv
target_data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20241111_combined_WGBS/
samples=$(find ${add_data_path} -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort)  

cleaned_file=$(mktemp)  
cat "$previous_data_path_list" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  

max_jobs=6
current_jobs() {  
    jobs -rp | wc -l  
}  

for sample in ${samples[@]}
do
    mkdir -p ${target_data_path}${sample}
    previous_data_path=$(awk -v sample="$sample" -F, '$1 == sample {print $2}' "$cleaned_file") 
    echo ${sample} in ${previous_data_path} 

    f1=${previous_data_path}/${sample}*/*${sample}*_R1*.fastq.gz
    f2=${add_data_path}${sample}/*${sample}*_R1*.fastq.gz
    f3=${target_data_path}${sample}/${sample}_R1.fastq.gz
    echo 'zcat' ${f1} ${f2} '| gzip -> ' ${f3}
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    zcat ${f1} ${f2} | gzip -> ${f3} &

    f1=${previous_data_path}/${sample}*/*${sample}*_R2*.fastq.gz
    f2=${add_data_path}${sample}/*${sample}*_R2*.fastq.gz
    f3=${target_data_path}${sample}/${sample}_R2.fastq.gz
    echo 'zcat' ${f1} ${f2} '| gzip -> ' ${f3}
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    zcat ${f1} ${f2} | gzip -> ${f3} &
done

wait
echo merge done
