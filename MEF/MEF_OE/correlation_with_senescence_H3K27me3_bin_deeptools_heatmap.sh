data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF_OE/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/MEF_OE/
antibodys=(H3K27me3 H2AK119ub1)
conditions=(Bmi1 Cbx2 Cbx7)
mkdir -p ${result_path}matrix_bin_heatmap/
mkdir -p ${result_path}plot_bin_heatmap/
quadrants=(first second third fourth)
peaks=()


for condition in ${conditions[@]}
do
    files=()
    for quadrant in ${quadrants[@]}
    do
        peak=${data_path}H3K27me3/MEF_${condition}/bed/H3K27me3_correlation_with_senescence_10kb_bin_${quadrant}.bed
        peaks+=("$peak")
    done
    for antibody in ${antibodys[@]}
    do
        echo ${antibody}
        search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff_batch.csv
        cleaned_file=$(mktemp)  
        cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
        samples=$(awk -F',' -v c="MEF_Vector" -v a="$antibody" 'NR > 1 && ($1 == c) && ($2 == a) {print $3}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a sample_array <<<"$samples"  
        echo ${sample_array[@]}
        for sample in ${sample_array[@]}
        do
            file=$(ls ${data_path}${antibody}/MEF_Vector/bw/${sample}n.bw)
            files+=("$file")
        done

        samples=$(awk -F',' -v c="MEF_$condition" -v a="$antibody" 'NR > 1 && ($1 == c) && ($2 == a) {print $3}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a sample_array <<<"$samples"  
        echo ${sample_array[@]}
        for sample in ${sample_array[@]}
        do
            file=$(ls ${data_path}${antibody}/MEF_${condition}/bw/${sample}n.bw)
            files+=("$file")
        done
    done
            
    samples=(NTY462 NTY463 NTY464 NTY465)
    for sample in ${samples[@]}
    do
        file=$(ls /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF/H3K27me3/bw/${sample}*_bs1000.bw)
        files+=("$file")
    done
    samples=(NTY466 NTY467 NTY468 NTY469)
    for sample in ${samples[@]}
    do
        file=$(ls /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF/H3K9me3/bw/${sample}*_bs1000.bw)
        files+=("$file")
    done

    computeMatrix scale-regions -S ${files[@]} -R ${peaks[@]} \
        --beforeRegionStartLength 10000 --startLabel Start --endLabel End \
        --regionBodyLength 10000 \
        --afterRegionStartLength 10000 \
        --numberOfProcessors 10 \
        --skipZeros -o ${result_path}matrix_bin_heatmap/MEF_${condition}.mat.gz &
    
done

for condition in ${conditions[@]}
do
    # plotHeatmap -m ${result_path}matrix_bin_heatmap/MEF_${condition}.mat.gz \
    #     -out ${result_path}plot_bin_heatmap/MEF_${condition}.png \
    #     --colorList 'white,blue' \
    #     --startLabel Start --endLabel End \
    #     --regionsLabel First Second Third Fourth \
    #     --zMin 0 --zMax 4 \
    #     --samplesLabel "Vector H3K27me3" "Vector H3K27me3" "${condition} H3K27me3" "${condition} H3K27me3" "Vector H2AK119ub1" "Vector H2AK119ub1" "${condition} H2AK119ub1" "${condition} H2AK119ub1" "p2 H3K27me3" "p2 H3K27me3" "p10 H3K27me3" "p10 H3K27me3" "p2 H3K9me3" "p2 H3K9me3" "p10 H3K9me3" "p10 H3K9me3" \
    #     --whatToShow 'heatmap and colorbar' &
    plotHeatmap -m ${result_path}matrix_bin_heatmap/MEF_${condition}.mat.gz \
        -out ${result_path}plot_bin_heatmap/MEF_${condition}2.png \
        --colorList 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'white,red'  'white,red'  'white,red'  'white,red' \
        --zMin 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \
        --zMax 4 4 4 4 4 4 4 4 4 4 4 4 1 1 1 1 \
        --startLabel Start --endLabel End \
        --regionsLabel First Second Third Fourth \
        --samplesLabel "Vector H3K27me3" "Vector H3K27me3" "${condition} H3K27me3" "${condition} H3K27me3" "Vector H2AK119ub1" "Vector H2AK119ub1" "${condition} H2AK119ub1" "${condition} H2AK119ub1" "p2 H3K27me3" "p2 H3K27me3" "p10 H3K27me3" "p10 H3K27me3" "p2 H3K9me3" "p2 H3K9me3" "p10 H3K9me3" "p10 H3K9me3" \
        --whatToShow 'heatmap and colorbar' &
done