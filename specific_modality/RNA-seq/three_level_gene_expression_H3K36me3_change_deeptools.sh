# conditions=(Stable Up Down low medium high)
conditions=(low medium)
# tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus ileum mammarygland iWAT)
tissues=(bladder pancreas spleen testis tongue uterus)
for tissue in ${tissues[@]}
do
    antibody=H3K36me3
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/${tissue}/${antibody}/
    RNA_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/${tissue}/
    search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/combined_analysis_enhancer/CUTTag_search_table.csv
    result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/${tissue}/
    mkdir -p ${result_path}three_level_gene_expression/
    mkdir -p ${result_path}three_level_gene_expression/matrix
    mkdir -p ${result_path}three_level_gene_expression/plot
    cleaned_file=$(mktemp)  
    cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
    young_samples=$(awk -F',' -v t="$tissue" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == "young") && ($2 == a) {print $3}' "$cleaned_file") 
    old_samples=$(awk -F',' -v t="$tissue" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == "old") && ($2 == a) {print $3}' "$cleaned_file") 
    IFS=$'\n' read -r -d '' -a young_array < <(echo "$young_samples" && printf '\0') 
    IFS=$'\n' read -r -d '' -a old_array < <(echo "$old_samples" && printf '\0') 
    # young1_bw=${data_path}/bw/${young_array[0]}*nodup.bw
    # young2_bw=${data_path}/bw/${young_array[1]}*nodup.bwread
    # old1_bw=${data_path}/bw/${old_array[0]}*nodup.bw
    # old2_bw=${data_path}/bw/${old_array[1]}*nodup.bw
    young1_bw=${data_path}/bw/${young_array[0]}*bs1000.bw
    young2_bw=${data_path}/bw/${young_array[1]}*bs1000.bw
    old1_bw=${data_path}/bw/${old_array[0]}*bs1000.bw
    old2_bw=${data_path}/bw/${old_array[1]}*bs1000.bw
    for condition in ${conditions[@]}
    do
        bed=${RNA_path}bed/${condition}_expression_gene.bed
        computeMatrix scale-regions -S $young1_bw $young2_bw $old1_bw $old2_bw -R $bed \
            --beforeRegionStartLength 10000 --startLabel TSS --endLabel TES \
            --regionBodyLength 10000 \
            --afterRegionStartLength 10000 \
            --numberOfProcessors 10 \
            --skipZeros -o ${result_path}three_level_gene_expression/matrix/${antibody}_in_${condition}_gene_expression_level.mat.gz
        
        plotProfile -m ${result_path}three_level_gene_expression/matrix/${antibody}_in_${condition}_gene_expression_level.mat.gz \
            --plotTitle "${tissue} ${antibody} in ${condition} expression gene" \
            --samplesLabel "Young1" "Young2" "Old1" "Old2" \
            --colors "#f38181" "#ff2e63" "#112d4e" "#3f72af" \
            --plotHeight 10 \
            --plotWidth 12 \
            --yMin 0 \
            --yMax 4 \
            --regionsLabel "Regions" \
            --yAxisLabel "Signal" \
            --legendLocation "upper-right" \
            --refPointLabel "Center" \
            --perGroup \
            -out ${result_path}three_level_gene_expression/plot/${antibody}_in_${condition}_gene_expression_level.pdf
    done
done