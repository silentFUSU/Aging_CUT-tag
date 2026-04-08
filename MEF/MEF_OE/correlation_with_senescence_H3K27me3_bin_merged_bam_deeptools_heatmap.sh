data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF_OE/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/MEF_OE/
# antibodys=(H3K27me3 H2AK119ub1)
antibodys=(H3K27me3)
# conditions=(Bmi1 Cbx2 Cbx7)
conditions=(Cbx7)
mkdir -p ${result_path}matrix_bin_heatmap/
mkdir -p ${result_path}plot_bin_heatmap/
# quadrants=(first second third fourth)
quadrants=(second fourth)

for condition in ${conditions[@]}
do
    files=()
    peaks=()
    for quadrant in ${quadrants[@]}
    do
        peak=${data_path}H3K27me3/MEF_${condition}/bed/H3K27me3_correlation_with_senescence_10kb_bin_${quadrant}.bed
        peaks+=("$peak")
    done
    # peak=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/Bmi1_GSE112460/peaks/Bmi1_peaks.bed
    # peaks+=("$peak")
    for antibody in ${antibodys[@]}
    do
        file=$(ls ${data_path}${antibody}/MEF_Vector/bw/${antibody}_MEF_Vector.bw)
        files+=("$file")
        file=$(ls ${data_path}${antibody}/MEF_${condition}/bw/${antibody}_MEF_${condition}.bw)
        files+=("$file")
    done
    file=$(ls /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF/H3K27me3/bw/young_bs1000.bw)
    files+=("$file")
    file=$(ls /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF/H3K27me3/bw/old_bs1000.bw)
    files+=("$file")

    # file=$(ls /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF/H3K9me3/bw/young_bs1000.bw)
    # files+=("$file")
    # file=$(ls /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF/H3K9me3/bw/old_bs1000.bw)
    # files+=("$file")
    
    file=$(ls /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/Bmi1_GSE112460/bigWig/Bmi1_logFE.sorted.nonneg.bw)
    files+=("$file")

    file=$(ls /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/Cbx7_GSE151901/bigWig/Cbx7_logFE.sorted.nonneg.bw)
    files+=("$file")

    computeMatrix scale-regions -S ${files[@]} -R ${peaks[@]} \
        --beforeRegionStartLength 50000 --startLabel Start --endLabel End \
        --regionBodyLength 10000 \
        --afterRegionStartLength 50000 \
        --missingDataAsZero \
        --numberOfProcessors 10 \
        --skipZeros -o ${result_path}matrix_bin_heatmap/MEF_${condition}_with_Bmi1_Cbx7_logFE_bw_50kb_II_IV.mat.gz &
    
    # computeMatrix scale-regions -S ${files[@]} -R ${peaks[@]} \
    #     --beforeRegionStartLength 10000 --startLabel Start --endLabel End \
    #     --regionBodyLength 10000 \
    #     --afterRegionStartLength 10000 \
    #     --missingDataAsZero \
    #     --numberOfProcessors 10 \
    #     --skipZeros -o ${result_path}matrix_bin_heatmap/MEF_${condition}_with_Bmi1_Cbx7_logFE_bw.mat.gz &
    
done
wait
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
    # plotHeatmap -m ${result_path}matrix_bin_heatmap/MEF_${condition}_with_Bmi1_peaks.mat.gz \
    #     -out ${result_path}plot_bin_heatmap/MEF_${condition}_with_Bmi1_peaks.pdf \
    #     --colorList 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'yellow,red'  'yellow,red' \
    #     --zMin 0 0 0 0 0 0 0 0  \
    #     --zMax 4 4 4 4 4 4 1 1 \
    #     --startLabel Start --endLabel End \
    #     --regionsLabel First Second Third Fourth "Bmi1 peaks" \
    #     --samplesLabel "Vector H3K27me3" "${condition} H3K27me3" "Vector H2AK119ub1" "${condition} H2AK119ub1" "p2 H3K27me3" "p10 H3K27me3" "p2 H3K9me3" "p10 H3K9me3" \
    #     --whatToShow 'heatmap and colorbar' &
    # plotHeatmap -m ${result_path}matrix_bin_heatmap/MEF_${condition}_with_Bmi1_Cbx7_logFE_bw_50kb.mat.gz \
    #     -out ${result_path}plot_bin_heatmap/MEF_${condition}_with_Bmi1_Cbx7_logFE_bw_50kb.pdf \
    #     --colorList '#f7fbff,#08306b' '#f7fbff,#08306b' '#f7fbff,#08306b' '#f7fbff,#08306b' '#f7fbff,#08306b' '#f7fbff,#08306b' '#ffffcc,#800026'  '#ffffcc,#800026' '#fcfbfd,#5e3c99,#1f0033' '#fcfbfd,#5e3c99,#1f0033' \
    #     --zMin 0 0 0 0 0 0 0 0 0 0 \
    #     --zMax 3 3 3 3 3 3 1 1 0.05 0.05 \
    #     --startLabel Start --endLabel End \
    #     --regionsLabel First Second Third Fourth \
    #     --samplesLabel "Vector H3K27me3" "${condition} H3K27me3" "Vector H2AK119ub1" "${condition} H2AK119ub1" "p2 H3K27me3" "p10 H3K27me3" "p2 H3K9me3" "p10 H3K9me3" "Bmi1" "Cbx7" \
    #     --whatToShow 'heatmap and colorbar' &
    
    # plotHeatmap -m ${result_path}matrix_bin_heatmap/MEF_${condition}_with_Bmi1_Cbx7_logFE_bw.mat.gz \
    #     -out ${result_path}plot_bin_heatmap/MEF_${condition}_with_Bmi1_Cbx7_logFE_bw.pdf \
    #     --colorList 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'white,blue' 'yellow,red'  'yellow,red' 'white,purple' 'white,purple'\
    #     --zMin 0 0 0 0 0 0 0 0 0 0 \
    #     --zMax 4 4 4 4 4 4 1 1 0.5 0.5 \
    #     --startLabel Start --endLabel End \
    #     --regionsLabel First Second Third Fourth \
    #     --samplesLabel "Vector H3K27me3" "${condition} H3K27me3" "Vector H2AK119ub1" "${condition} H2AK119ub1" "p2 H3K27me3" "p10 H3K27me3" "p2 H3K9me3" "p10 H3K9me3" "Bmi1" "Cbx7" \
    #     --whatToShow 'heatmap and colorbar' &
    plotHeatmap -m ${result_path}matrix_bin_heatmap/MEF_${condition}_with_Bmi1_Cbx7_logFE_bw_50kb_II_IV.mat.gz \
        -out ${result_path}plot_bin_heatmap/MEF_${condition}_with_Bmi1_Cbx7_logFE_bw_50kb_II_IV.pdf \
        --colorList '#f7fbff,#08306b' '#f7fbff,#08306b' '#f7fbff,#08306b' '#f7fbff,#08306b' '#fcfbfd,#5e3c99,#1f0033' '#fcfbfd,#5e3c99,#1f0033' \
        --zMin 0 0 0 0 0 0 0 \
        --zMax 3 3 3 3 0.05 0.05 \
        --yMin 0 0 0 0 -0.01 -0.01 \
        --yMax 3 3 3 3 0.04 0.04 \
        --startLabel Start --endLabel End \
        --regionsLabel Second Fourth \
        --samplesLabel "Vector H3K27me3" "${condition} H3K27me3"  "p2 H3K27me3" "p10 H3K27me3"  "Bmi1" "Cbx7"  &
    # plotProfile -m ${result_path}matrix_bin_heatmap/MEF_${condition}_with_Bmi1_Cbx7_logFE_bw_50kb_II_IV.mat.gz \
    #         --plotTitle "${tissue} ${antibody}" \
    #         --samplesLabel "Young"  "Old" \
    #         --colors "#e64b35" "#3c5488" \
    #         --plotHeight 10 \
    #         --plotWidth 10 \
    #         --regionsLabel "Regions" \
    #         --yAxisLabel "Signal" \
    #         --legendLocation "upper-right" \
    #         --refPointLabel "Center" \
    #         --perGroup \
    #         --startLabel Start --endLabel End \
    #         --yMax 0.5 \
    #         -out ${result_path}plot_merge/${tissue}_${antibody}_change_in_H3K27me3_domain.pdf
done