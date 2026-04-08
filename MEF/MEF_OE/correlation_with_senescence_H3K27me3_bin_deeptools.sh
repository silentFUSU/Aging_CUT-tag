data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF_OE/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/MEF_OE/
antibodys=(H3K27me3 H2AK119ub1)
conditions=(MEF_Bmi1 MEF_Cbx2 MEF_Cbx7)
mkdir -p ${result_path}matrix_bin/
mkdir -p ${result_path}plot_bin/
quadrants=(first second third fourth)
for antibody in ${antibodys[@]}
do
    for condition in ${conditions[@]}
    do
        search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff_batch.csv
        cleaned_file=$(mktemp)  
        cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
        samples=$(awk -F',' -v c="MEF_Vector" -v a="$antibody" 'NR > 1 && ($1 == c) && ($2 == a) {print $3}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a sample_array <<<"$samples"  
        echo ${sample_array[@]}
        files=()
        for sample in ${sample_array[@]}
        do
            file=$(ls ${data_path}${antibody}/MEF_Vector/bw/${sample}n.bw)
            files+=("$file")
        done

        samples=$(awk -F',' -v c="$condition" -v a="$antibody" 'NR > 1 && ($1 == c) && ($2 == a) {print $3}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a sample_array <<<"$samples"  
        echo ${sample_array[@]}
        for sample in ${sample_array[@]}
        do
            file=$(ls ${data_path}${antibody}/${condition}/bw/${sample}n.bw)
            files+=("$file")
        done

        for quadrant in ${quadrants[@]}
        do

            peak=${data_path}H3K27me3/${condition}/bed/H3K27me3_correlation_with_senescence_10kb_bin_${quadrant}.bed
    
            # computeMatrix scale-regions -S ${files[@]} -R $peak \
            #     --beforeRegionStartLength 10000 --startLabel Start --endLabel End \
            #     --regionBodyLength 10000 \
            #     --afterRegionStartLength 10000 \
            #     --numberOfProcessors 10 \
            #     --skipZeros -o ${result_path}matrix_bin/${antibody}_${condition}_${quadrant}.mat.gz
            
            plotProfile -m ${result_path}matrix_bin/${antibody}_${condition}_${quadrant}.mat.gz \
                --plotTitle "${antibody} ${condition} ${quadrant}" \
                --colors "#f38181" "#ff2e63" "#112d4e" "#3f72af" \
                --plotHeight 10 \
                --plotWidth 10 \
                --regionsLabel "Regions" \
                --yAxisLabel "Signal" \
                --legendLocation "upper-right" \
                --refPointLabel "Center" \
                --perGroup \
                --startLabel Start --endLabel End \
                --yMin 0.1 --yMax 2.5 \
                -out ${result_path}plot_bin/${antibody}_${condition}_${quadrant}.pdf 
        done
    done
done

for condition in ${conditions[@]}
do
    files=()
    samples=(NTY462 NTY463)

    for sample in ${samples[@]}
    do
        file=$(ls /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF/H3K27me3/bw/${sample}*_bs1000.bw)
        files+=("$file")
    done
    samples=(NTY464 NTY465)
    for sample in ${samples[@]}
    do
        file=$(ls /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF/H3K27me3/bw/${sample}*_bs1000.bw)
        files+=("$file")
    done

    for quadrant in ${quadrants[@]}
    do
        peak=${data_path}H3K27me3/${condition}/bed/H3K27me3_correlation_with_senescence_10kb_bin_${quadrant}.bed
        # computeMatrix scale-regions -S ${files[@]} -R $peak \
        #     --beforeRegionStartLength 10000 --startLabel Start --endLabel End \
        #     --regionBodyLength 10000 \
        #     --afterRegionStartLength 10000 \
        #     --numberOfProcessors 10 \
        #     --skipZeros -o ${result_path}matrix_bin/senescence_H3K27me3_in_${condition}_${quadrant}.mat.gz 
        
        plotProfile -m ${result_path}matrix_bin/senescence_H3K27me3_in_${condition}_${quadrant}.mat.gz  \
            --plotTitle "senescence H3K27me3 in ${condition} ${quadrant}" \
            --colors "#f38181" "#ff2e63" "#112d4e" "#3f72af" \
            --plotHeight 10 \
            --plotWidth 10 \
            --regionsLabel "Regions" \
            --yAxisLabel "Signal" \
            --legendLocation "upper-right" \
            --refPointLabel "Center" \
            --perGroup \
            --startLabel Start --endLabel End \
            --yMin 0.1 --yMax 2.5 \
            -out ${result_path}plot_bin/senescence_H3K27me3_in_${condition}_${quadrant}.pdf 
    done
done