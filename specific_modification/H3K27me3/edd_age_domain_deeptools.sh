tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum bonemarrow ileum heart thymus stomach skin aorta tongue bladder CB jejunum uterus ovary BAT iWAT mammarygland)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/H3K27me3_domain/
mkdir -p $result_path
max_jobs=4
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do
    search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff_batch.csv
    mkdir -p ${result_path}matrix_median
    mkdir -p ${result_path}plot
    cleaned_file=$(mktemp)  
    cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
    young_samples=$(awk -F',' -v t="$tissue" -v a="H3K27me3" 'NR > 1 && ($1 == t) && ($5 == "3m") && ($2 == a) {print $3}' "$cleaned_file") 
    old_samples=$(awk -F',' -v t="$tissue" -v a="H3K27me3" 'NR > 1 && ($1 == t) && ($5 == "24m") && ($2 == a) {print $3}' "$cleaned_file") 
    IFS=$'\n' read -r -d '' -a young_array < <(echo "$young_samples" && printf '\0') 
    IFS=$'\n' read -r -d '' -a old_array < <(echo "$old_samples" && printf '\0') 
    young1_bw=${data_path}${tissue}/H3K27me3/bw/${young_array[0]}*nodup.bw
    young2_bw=${data_path}${tissue}/H3K27me3/bw/${young_array[1]}*nodup.bw
    old1_bw=${data_path}${tissue}/H3K27me3/bw/${old_array[0]}*nodup.bw
    old2_bw=${data_path}${tissue}/H3K27me3/bw/${old_array[1]}*nodup.bw
    bed=${data_path}${tissue}/H3K27me3/peaks/edd/edd_peaks_fdr05.bed
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    computeMatrix scale-regions -S $young1_bw $young1_bw $old1_bw $old2_bw -R $bed \
        --beforeRegionStartLength 10000 --startLabel Start --endLabel End \
        --regionBodyLength 10000 \
        --afterRegionStartLength 10000 \
        --numberOfProcessors 10 \
        --averageTypeBins median \
        --skipZeros -o ${result_path}matrix_median/${tissue}_H3K27me3_change_in_H3K27me3_domain.mat.gz &
    
    search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff_batch.csv
    mkdir -p ${result_path}matrix_median
    mkdir -p ${result_path}plot
    cleaned_file=$(mktemp)  
    cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
    young_samples=$(awk -F',' -v t="$tissue" -v a="H3K9me3" 'NR > 1 && ($1 == t) && ($5 == "3m") && ($2 == a) {print $3}' "$cleaned_file") 
    old_samples=$(awk -F',' -v t="$tissue" -v a="H3K9me3" 'NR > 1 && ($1 == t) && ($5 == "24m") && ($2 == a) {print $3}' "$cleaned_file") 
    IFS=$'\n' read -r -d '' -a young_array < <(echo "$young_samples" && printf '\0') 
    IFS=$'\n' read -r -d '' -a old_array < <(echo "$old_samples" && printf '\0') 
    young1_bw=${data_path}${tissue}/H3K9me3/bw/${young_array[0]}*nodup.bw
    young2_bw=${data_path}${tissue}/H3K9me3/bw/${young_array[1]}*nodup.bw
    old1_bw=${data_path}${tissue}/H3K9me3/bw/${old_array[0]}*nodup.bw
    old2_bw=${data_path}${tissue}/H3K9me3/bw/${old_array[1]}*nodup.bw
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    computeMatrix scale-regions -S $young1_bw $young1_bw $old1_bw $old2_bw -R $bed \
        --beforeRegionStartLength 10000 --startLabel Start --endLabel End \
        --regionBodyLength 10000 \
        --afterRegionStartLength 10000 \
        --numberOfProcessors 10 \
        --averageTypeBins median \
        --skipZeros -o ${result_path}matrix_median/${tissue}_H3K9me3_change_in_H3K27me3_domain.mat.gz &
    
done