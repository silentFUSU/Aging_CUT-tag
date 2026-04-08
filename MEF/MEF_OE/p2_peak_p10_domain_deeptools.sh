data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF_OE/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/MEF_OE/
antibodys=(H3K27me3 H2AK119ub1)
conditions=(MEF_Bmi1 MEF_Cbx2 MEF_Cbx7)
window_size=5000
gap_size=10000
e_value=100
mkdir -p ${result_path}matrix/
mkdir -p ${result_path}plot/

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
        peak=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed
        
        computeMatrix scale-regions -S ${files[@]} -R $peak \
            --beforeRegionStartLength 10000 --startLabel Start --endLabel End \
            --regionBodyLength 10000 \
            --afterRegionStartLength 10000 \
            --numberOfProcessors 10 \
            --skipZeros -o ${result_path}matrix/${antibody}_${condition}_p2_peaks.mat.gz
        
        plotProfile -m ${result_path}matrix/${antibody}_${condition}_p2_peaks.mat.gz \
            --plotTitle "${antibody} ${condition}" \
            --colors "#f38181" "#ff2e63" "#112d4e" "#3f72af" \
            --plotHeight 10 \
            --plotWidth 10 \
            --regionsLabel "Regions" \
            --yAxisLabel "Signal" \
            --legendLocation "upper-right" \
            --refPointLabel "Center" \
            --perGroup \
            --startLabel Start --endLabel End \
            -out ${result_path}plot/${antibody}_${condition}_p2_peaks.pdf 

        domain=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF/H3K27me3/peaks/edd/edd_peaks_fdr05.bed
        computeMatrix scale-regions -S ${files[@]} -R $domain \
            --beforeRegionStartLength 10000 --startLabel Start --endLabel End \
            --regionBodyLength 10000 \
            --afterRegionStartLength 10000 \
            --numberOfProcessors 10 \
            --skipZeros -o ${result_path}matrix/${antibody}_${condition}_p10_domain.mat.gz
        plotProfile -m ${result_path}matrix/${antibody}_${condition}_p10_domain.mat.gz \
            --plotTitle "${antibody} ${condition}" \
            --colors "#f38181" "#ff2e63" "#112d4e" "#3f72af" \
            --plotHeight 10 \
            --plotWidth 10 \
            --regionsLabel "Regions" \
            --yAxisLabel "Signal" \
            --legendLocation "upper-right" \
            --refPointLabel "Center" \
            --perGroup \
            --startLabel Start --endLabel End \
            -out ${result_path}plot/${antibody}_${condition}_p10_domain.pdf 
    done
done