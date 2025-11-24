tissue=mammarygland
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/H3K9me3/recursion_peaks/
mkdir -p ${result_path}TE_with_histone/
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
antibody=H3K27me3
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff.csv

kmeans=(kmeans1 kmeans2 kmeans3 kmeans4)
for tissue in ${tissues[@]}
do
    for kmean in ${kmeans[@]}
    do
        bed=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/${kmean}_uinon_recursion_peaks.bed
        mkdir -p ${result_path}${kmean}
        mkdir -p ${result_path}${kmean}/matrix
        mkdir -p ${result_path}${kmean}/plot
        search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff_batch.csv
        cleaned_file=$(mktemp)  
        cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
        young_samples_H3K27me3_rep1=$(awk -F',' -v t="$tissue" -v a="$antibody" -v b="batch1" 'NR > 1 && ($1 == t) && ($5 == "3m") && ($2 == a) && ($6 == b) {print $3}' "$cleaned_file") 
        old_samples_H3K27me3_rep1=$(awk -F',' -v t="$tissue" -v a="$antibody" -v b="batch1" 'NR > 1 && ($1 == t) && ($5 == "24m") && ($2 == a) && ($6 == b) {print $3}' "$cleaned_file") 
        young_samples_H3K27me3_rep2=$(awk -F',' -v t="$tissue" -v a="$antibody" -v b="batch2" 'NR > 1 && ($1 == t) && ($5 == "3m") && ($2 == a) && ($6 == b) {print $3}' "$cleaned_file") 
        old_samples_H3K27me3_rep2=$(awk -F',' -v t="$tissue" -v a="$antibody" -v b="batch2" 'NR > 1 && ($1 == t) && ($5 == "24m") && ($2 == a) && ($6 == b) {print $3}' "$cleaned_file")
        
        young_samples_H3K9me3_rep1=$(awk -F',' -v t="$tissue" -v a="H3K9me3" -v b="batch1" 'NR > 1 && ($1 == t) && ($5 == "3m") && ($2 == a) && ($6 == b) {print $3}' "$cleaned_file") 
        old_samples_H3K9me3_rep1=$(awk -F',' -v t="$tissue" -v a="H3K9me3" -v b="batch1" 'NR > 1 && ($1 == t) && ($5 == "24m") && ($2 == a) && ($6 == b) {print $3}' "$cleaned_file") 
        young_samples_H3K9me3_rep2=$(awk -F',' -v t="$tissue" -v a="H3K9me3" -v b="batch2" 'NR > 1 && ($1 == t) && ($5 == "3m") && ($2 == a) && ($6 == b) {print $3}' "$cleaned_file") 
        old_samples_H3K9me3_rep2=$(awk -F',' -v t="$tissue" -v a="H3K9me3" -v b="batch2" 'NR > 1 && ($1 == t) && ($5 == "24m") && ($2 == a) && ($6 == b) {print $3}' "$cleaned_file")
        
        young1_H3K27me3_bw=${data_path}${tissue}/${antibody}/bw/${young_samples_H3K27me3_rep1}*nodup.bw
        young2_H3K27me3_bw=${data_path}${tissue}/${antibody}/bw/${young_samples_H3K27me3_rep2}*nodup.bw
        old1_H3K27me3_bw=${data_path}${tissue}/${antibody}/bw/${old_samples_H3K27me3_rep1}*nodup.bw
        old2_H3K27me3_bw=${data_path}${tissue}/${antibody}/bw/${old_samples_H3K27me3_rep2}*nodup.bw

        young1_H3K9me3_bw=${data_path}${tissue}/H3K9me3/bw/${young_samples_H3K9me3_rep1}*nodup.bw
        young2_H3K9me3_bw=${data_path}${tissue}/H3K9me3/bw/${young_samples_H3K9me3_rep2}*nodup.bw
        old1_H3K9me3_bw=${data_path}${tissue}/H3K9me3/bw/${old_samples_H3K9me3_rep1}*nodup.bw
        old2_H3K9me3_bw=${data_path}${tissue}/H3K9me3/bw/${old_samples_H3K9me3_rep2}*nodup.bw

    
        computeMatrix scale-regions -S $young1_H3K27me3_bw $young2_H3K27me3_bw  $young1_H3K9me3_bw $young2_H3K9me3_bw -R $bed \
            --beforeRegionStartLength 1000 --startLabel Start --endLabel End \
            --regionBodyLength 1000 \
            --afterRegionStartLength 1000 \
            --numberOfProcessors 10 \
            --skipZeros -o ${result_path}${kmean}/matrix/${tissue}_${kmean}_${antibody}_H3K9me3.mat.gz

        plotProfile -m ${result_path}${kmean}/matrix/${tissue}_${kmean}_${antibody}_H3K9me3.mat.gz \
            --plotTitle "${tissue} ${antibody} H3K9me3 in ${kmean} regions" \
            --samplesLabel "Young1_${antibody}" "Young2_${antibody}" "Young1_H3K9me3" "Young2_H3K9me3" \
            --colors "#f38181" "#ff2e63" "#112d4e" "#3f72af" \
            --plotHeight 10 \
            --plotWidth 12 \
            --regionsLabel "Regions" \
            --yAxisLabel "Signal" \
            --legendLocation "upper-right" \
            --refPointLabel "Center" \
            --perGroup \
            -out  ${result_path}${kmean}/plot/${tissue}_${kmean}_${antibody}_H3K9me3.pdf
    done
done