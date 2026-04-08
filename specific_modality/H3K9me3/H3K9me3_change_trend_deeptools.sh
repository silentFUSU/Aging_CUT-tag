tissues=(testis tongue)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
antibody=H3K9me3
for tissue in ${tissues[@]}
do
    result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/${tissue}/H3K9me3/
    mkdir ${result_path}
    mkdir -p ${result_path}matrix
    mkdir -p ${result_path}plot
    search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff.csv
    cleaned_file=$(mktemp)  
    cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
    young_samples=$(awk -F',' -v t="$tissue" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == "3m") && ($2 == a) {print $3}' "$cleaned_file") 
    old_samples=$(awk -F',' -v t="$tissue" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == "24m") && ($2 == a) {print $3}' "$cleaned_file") 
    IFS=$'\n' read -r -d '' -a young_array < <(echo "$young_samples" && printf '\0') 
    IFS=$'\n' read -r -d '' -a old_array < <(echo "$old_samples" && printf '\0') 
    young1_bw=${data_path}${tissue}/${antibody}/bw/${young_array[0]}*bs1000.bw
    young2_bw=${data_path}${tissue}/${antibody}/bw/${young_array[1]}*bs1000.bw
    old1_bw=${data_path}${tissue}/${antibody}/bw/${old_array[0]}*bs1000.bw
    old2_bw=${data_path}${tissue}/${antibody}/bw/${old_array[1]}*bs1000.bw
    bed=${data_path}all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_bedtools_filtered.bed
    
    computeMatrix scale-regions -S $young1_bw $young2_bw $old1_bw $old2_bw -R $bed \
        --beforeRegionStartLength 10000 --startLabel Start --endLabel End \
        --regionBodyLength 10000 \
        --afterRegionStartLength 10000 \
        --numberOfProcessors 10 \
        --skipZeros -o ${result_path}matrix/${tissue}_H3K9me3_in_H3K9me3-W5000-G10000-E100_bedtools_filter_peaks.mat.gz &
    
    plotProfile -m ${result_path}matrix/${tissue}_H3K9me3_in_H3K9me3-W5000-G10000-E100_bedtools_filter_peaks.mat.gz \
        --plotTitle "${tissue} H3K9me3" \
        --samplesLabel "Young1" "Young2" "Old1" "Old2" \
        --colors "#f38181" "#ff2e63" "#112d4e" "#3f72af" \
        --plotHeight 10 \
        --plotWidth 10 \
        --regionsLabel "Regions" \
        --yAxisLabel "Signal" \
        --legendLocation "upper-right" \
        --refPointLabel "Center" \
        --perGroup \
        --startLabel Start --endLabel End \
        -out ${result_path}plot/${tissue}_H3K9me3_in_H3K9me3-W5000-G10000-E100_bedtools_filter_peaks.pdf &
done