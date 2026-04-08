data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/diff/H3K9me3/
mkdir ${result_path}recursion_peaks_kmeans
mkdir ${result_path}recursion_peaks_kmeans/matrix
mkdir ${result_path}recursion_peaks_kmeans/plot
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
antibodys=(H3K27me3)
kmeans=(kmeans1 kmeans2 kmeans3 kmeans4)
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff.csv

for kmean in ${kmeans[@]}
do
    bed=${data_path}all/H3K9me3/recursion_peaks_diff_table/bed/${kmean}_uinon_recursion_peaks.bed
    for antibody in ${antibodys[@]}
    do
        for tissue in ${tissues[@]}
        do
        search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff_batch.csv
        cleaned_file=$(mktemp)  
        cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
        young_samples=$(awk -F',' -v t="$tissue" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == "3m") && ($2 == a) {print $3}' "$cleaned_file") 
        old_samples=$(awk -F',' -v t="$tissue" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == "24m") && ($2 == a) {print $3}' "$cleaned_file") 
        IFS=$'\n' read -r -d '' -a young_array < <(echo "$young_samples" && printf '\0') 
        IFS=$'\n' read -r -d '' -a old_array < <(echo "$old_samples" && printf '\0') 
        young1_bw=${data_path}${tissue}/${antibody}/bw/${young_array[0]}*nodup.bw
        young2_bw=${data_path}${tissue}/${antibody}/bw/${young_array[1]}*nodup.bw
        old1_bw=${data_path}${tissue}/${antibody}/bw/${old_array[0]}*nodup.bw
        old2_bw=${data_path}${tissue}/${antibody}/bw/${old_array[1]}*nodup.bw
        # young1_bw=${data_path}${tissue}/${antibody}/bw/${young_array[0]}*bs1000.bw
        # young2_bw=${data_path}${tissue}/${antibody}/bw/${young_array[1]}*bs1000.bw
        # old1_bw=${data_path}${tissue}/${antibody}/bw/${old_array[0]}*bs1000.bw
        # old2_bw=${data_path}${tissue}/${antibody}/bw/${old_array[1]}*bs1000.bw
        computeMatrix scale-regions -S $young1_bw $young2_bw $old1_bw $old2_bw -R $bed \
            --beforeRegionStartLength 1000 --startLabel Start --endLabel End \
            --regionBodyLength 1000 \
            --afterRegionStartLength 1000 \
            --numberOfProcessors 10 \
            --skipZeros -o ${result_path}recursion_peaks_kmeans/matrix/${tissue}_${antibody}_change_in_H3K9me3_recursion_peaks_${kmean}_regions.mat.gz

        plotProfile -m ${result_path}recursion_peaks_kmeans/matrix/${tissue}_${antibody}_change_in_H3K9me3_recursion_peaks_${kmean}_regions.mat.gz \
            --plotTitle "${tissue} ${antibody} in H3K9me3 ${kmean} regions" \
            --samplesLabel "Young1" "Young2" "Old1" "Old2" \
            --colors "#f38181" "#ff2e63" "#112d4e" "#3f72af" \
            --plotHeight 10 \
            --plotWidth 12 \
            --regionsLabel "Regions" \
            --yAxisLabel "Signal" \
            --legendLocation "upper-right" \
            --refPointLabel "Center" \
            --perGroup \
            -out ${result_path}recursion_peaks_kmeans/plot/${tissue}_${antibody}_change_in_H3K9me3_recursion_peaks_${kmean}_regions.pdf
        done
    done
done