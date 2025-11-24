# tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum bonemarrow ileum heart thymus stomach skin aorta tongue bladder CB jejunum uterus ovary BAT iWAT mammarygland)
tissues=(cecum colon jejunum ileum)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/H3K27me3_H3K9me3/H3K27me3_change_in_H3K9me3_peaks/
for tissue in ${tissues[@]}
do
    search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff.csv
    mkdir -p ${result_path}matrix
    mkdir -p ${result_path}plot
    cleaned_file=$(mktemp)  
    cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
    young_samples=$(awk -F',' -v t="$tissue" -v a="H3K27me3" 'NR > 1 && ($1 == t) && ($5 == "3m") && ($2 == a) {print $3}' "$cleaned_file") 
    old_samples=$(awk -F',' -v t="$tissue" -v a="H3K27me3" 'NR > 1 && ($1 == t) && ($5 == "24m") && ($2 == a) {print $3}' "$cleaned_file") 
    IFS=$'\n' read -r -d '' -a young_array < <(echo "$young_samples" && printf '\0') 
    IFS=$'\n' read -r -d '' -a old_array < <(echo "$old_samples" && printf '\0') 
    young1_bw=${data_path}${tissue}/H3K27me3/bw/${young_array[0]}*bs1000.bw
    young2_bw=${data_path}${tissue}/H3K27me3/bw/${young_array[1]}*bs1000.bw
    old1_bw=${data_path}${tissue}/H3K27me3/bw/${old_array[0]}*bs1000.bw
    old2_bw=${data_path}${tissue}/H3K27me3/bw/${old_array[1]}*bs1000.bw
    H3K27me3_peaks=${data_path}${tissue}/H3K27me3/bed/H3K27me3_young_merge-W1000-G3000-E100_compress.bed
    H3K9me3_peaks=${data_path}${tissue}/H3K9me3/bed/H3K9me3_young_merge-W1000-G3000-E100_compress.bed
    bedtools intersect -a ${H3K9me3_peaks} -b ${H3K27me3_peaks}  | cut -f 1-3 > ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3_H3K27me3_peak_overlap_region.bed
    /usr/local/lib64/R/bin/Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/H3K27me3_H3K9me3_relationship/Get_H3K9me3_peak_overlap_or_not_H3K27me3_peak.R $tissue
        
    bedtools merge -i ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3_not_overlap_with_H3K27me3_peaks.bed -d 30000 | \
    awk 'BEGIN{OFS="\t"} {print $0, "peak" NR}' > ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3_not_overlap_with_H3K27me3_peaks_compress.bed
    bed=${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3_not_overlap_with_H3K27me3_peaks_compress.bed
    computeMatrix scale-regions -S $young1_bw $young2_bw $old1_bw $old2_bw -R $bed \
        --beforeRegionStartLength 1000000 --startLabel Start --endLabel End \
        --regionBodyLength 50000 \
        --afterRegionStartLength 1000000 \
        --numberOfProcessors 10 \
        --skipZeros -o ${result_path}matrix/${tissue}_H3K27me3_change_in_H3K9me3_peaks_not_overlap_with_H3K27me3_peaks.mat.gz &
    
    bedtools merge -i ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3_overlap_with_H3K27me3_peaks.bed -d 30000 | \
    awk 'BEGIN{OFS="\t"} {print $0, "peak" NR}' > ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3_overlap_with_H3K27me3_peaks_compress.bed
    bed=${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3_overlap_with_H3K27me3_peaks_compress.bed
    computeMatrix scale-regions -S $young1_bw $young2_bw $old1_bw $old2_bw -R $bed \
        --beforeRegionStartLength 1000000 --startLabel Start --endLabel End \
        --regionBodyLength 50000 \
        --afterRegionStartLength 1000000 \
        --numberOfProcessors 10 \
        --skipZeros -o ${result_path}matrix/${tissue}_H3K27me3_change_in_H3K9me3_peaks_overlap_with_H3K27me3_peaks.mat.gz &
    
    wait

    plotProfile -m ${result_path}matrix/${tissue}_H3K27me3_change_in_H3K9me3_peaks_not_overlap_with_H3K27me3_peaks.mat.gz \
        --plotTitle "${tissue} H3K27me3 change in H3K9me3 peaks" \
        --samplesLabel "Young1" "Young2" "Old1" "Old2" \
        --colors "#f38181" "#ff2e63" "#112d4e" "#3f72af" \
        --plotHeight 10 \
        --plotWidth 30 \
        --regionsLabel "Regions" \
        --yAxisLabel "Signal" \
        --legendLocation "upper-right" \
        --refPointLabel "Center" \
        --perGroup \
        --startLabel Start --endLabel End \
        -out ${result_path}plot/${tissue}_H3K27me3_change_in_H3K9me3_peaks_not_overlap_with_H3K27me3_peaks.pdf &
    
    plotProfile -m ${result_path}matrix/${tissue}_H3K27me3_change_in_H3K9me3_peaks_overlap_with_H3K27me3_peaks.mat.gz \
        --plotTitle "${tissue} H3K27me3 change in H3K9me3 peaks" \
        --samplesLabel "Young1" "Young2" "Old1" "Old2" \
        --colors "#f38181" "#ff2e63" "#112d4e" "#3f72af" \
        --plotHeight 10 \
        --plotWidth 30 \
        --regionsLabel "Regions" \
        --yAxisLabel "Signal" \
        --legendLocation "upper-right" \
        --refPointLabel "Center" \
        --perGroup \
        --startLabel Start --endLabel End \
        -out ${result_path}plot/${tissue}_H3K27me3_change_in_H3K9me3_peaks_overlap_with_H3K27me3_peaks.pdf &
done
