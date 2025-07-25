class=ERVK_kmeans1
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/RNA/TE/
mkdir -p ${result_path}TE_with_histone/
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
TE_region=/storage/zhangyanxiaoLab/suzhuojie/ref_data/TE_reference/mm10_${class}_regions.bed
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
antibodys=(H3K9me3 H3K27me3 H3K36me3 H3K27ac H3K4me1 H3K4me3)
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff.csv
for antibody in ${antibodys[@]}
do
    mkdir -p ${result_path}TE_with_histone/${antibody}/
    mkdir -p ${result_path}TE_with_histone/${antibody}/matrix
    mkdir -p ${result_path}TE_with_histone/${antibody}/plot
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
        computeMatrix scale-regions -S $young1_bw $young2_bw $old1_bw $old2_bw -R $TE_region \
            --beforeRegionStartLength 1000 --startLabel Start --endLabel End \
            --regionBodyLength 1000 \
            --afterRegionStartLength 1000 \
            --numberOfProcessors 10 \
            --skipZeros -o ${result_path}TE_with_histone/${antibody}/matrix/${tissue}_${antibody}_change_in_${class}_regions.mat.gz
        
        plotProfile -m ${result_path}TE_with_histone/${antibody}/matrix/${tissue}_${antibody}_change_in_${class}_regions.mat.gz \
            --plotTitle "${tissue} ${antibody} in ${class} regions" \
            --samplesLabel "Young1" "Young2" "Old1" "Old2" \
            --colors "#f38181" "#ff2e63" "#112d4e" "#3f72af" \
            --plotHeight 10 \
            --plotWidth 12 \
            --regionsLabel "Regions" \
            --yAxisLabel "Signal" \
            --legendLocation "upper-right" \
            --refPointLabel "Center" \
            --perGroup \
            -out ${result_path}TE_with_histone/${antibody}/plot/${tissue}_${antibody}_change_in_${class}_regions.pdf
    done
done