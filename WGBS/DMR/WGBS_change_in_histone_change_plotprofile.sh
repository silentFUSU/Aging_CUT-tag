tissue=$1
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/${tissue}/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/WGBS/${tissue}/
CUTTag_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/WGBS_search_table.csv
mkdir -p ${result_path}WGBS_change_in_histone_change
mkdir -p ${result_path}WGBS_change_in_histone_change/matrix
mkdir -p ${result_path}WGBS_change_in_histone_change/plot
cleaned_file=$(mktemp)  
cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  

young_samples=$(awk -F',' -v t="$tissue"  'NR > 1 && ($1 == t) && ($5 == "3M") {print $3}' "$cleaned_file") 
old_samples=$(awk -F',' -v t="$tissue"  'NR > 1 && ($1 == t) && ($5 == "24M") {print $3}' "$cleaned_file") 
IFS=$'\n' read -r -d '' -a young_array < <(echo "$young_samples" && printf '\0') 
IFS=$'\n' read -r -d '' -a old_array < <(echo "$old_samples" && printf '\0') 
young1_bw=${data_path}/bw/${young_array[0]}*CpG.bw
young2_bw=${data_path}/bw/${young_array[1]}*CpG.bw
old1_bw=${data_path}/bw/${old_array[0]}*CpG.bw
old2_bw=${data_path}/bw/${old_array[1]}*CpG.bw

antibodys=(H3K27me3 H3K9me3 H3K36me3 H3K27ac H3K4me3 H3K4me1)
conditions=(up down)

for antibody in ${antibodys[@]}
do
    for condition in ${conditions[@]}
    do
        if [[ "$antibody" == "H3K27ac" || "$antibody" == "H3K4me3" || "$antibody" == "H3K4me1" ]]; then  
            bin_size="1kb"  
        else  
            bin_size="10kb"  
        fi  
        bed=${CUTTag_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_bins_diff_after_remove_batch_effect_${condition}.bed
        computeMatrix scale-regions -S $young1_bw $young2_bw $old1_bw $old2_bw -R $bed \
            --beforeRegionStartLength 10000 --startLabel start --endLabel end \
            --regionBodyLength 10000 \
            --afterRegionStartLength 10000 \
            --numberOfProcessors 10 \
            --skipZeros -o ${result_path}WGBS_change_in_histone_change/matrix/WGBS_change_in_${antibody}_${condition}.mat.gz

        plotProfile -m ${result_path}WGBS_change_in_histone_change/matrix/WGBS_change_in_${antibody}_${condition}.mat.gz \
            --plotTitle "WGBS Change in ${antibody} ${condition}" \
            --samplesLabel "Young1" "Young2" "Old1" "Old2" \
            --colors "blue" "blue" "red" "red" \
            --plotHeight 10 \
            --plotWidth 12 \
            --regionsLabel "Regions" \
            --yAxisLabel "Signal" \
            --legendLocation "upper-right" \
            --refPointLabel "Center" \
            --perGroup \
            -out ${result_path}WGBS_change_in_histone_change/plot/WGBS_change_in_${antibody}_${condition}.pdf
        /usr/local/lib64/R/bin/Rscript ~/projects/Aging_CUT_Tag/code/WGBS/DMR/WGBS_change_in_histone_change_plotprofile.R ${tissue} ${antibody} ${condition} &
    done
done