data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data_transposon2/20250908_MWY_CUTTAG/
kmeans=(kmeans1 kmeans2 kmeans3 kmeans4)
 for kmean in ${kmeans[@]}
    do
        bed=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/${kmean}_uinon_recursion_peaks.bed
        M2_bw=${data_path}bigWig/MWY923-1*nodup.bw
        M3_bw=${data_path}bigWig/MWY923-2*nodup.bw
        M10_bw=${data_path}bigWig/MWY923-3*nodup.bw
        M15_bw=${data_path}bigWig/MWY923-4*nodup.bw
        M19_bw=${data_path}bigWig/MWY923-5*nodup.bw
        M24_bw=${data_path}bigWig/MWY923-6*nodup.bw
        computeMatrix scale-regions -S $M2_bw $M3_bw $M10_bw $M15_bw $M19_bw $M24_bw -R $bed \
            --beforeRegionStartLength 1000 --startLabel Start --endLabel End \
            --regionBodyLength 1000 \
            --afterRegionStartLength 1000 \
            --numberOfProcessors 10 \
            --skipZeros -o ${data_path}matrix/${kmean}_H3K27me3.mat.gz
         plotProfile -m ${data_path}matrix/${kmean}_H3K27me3.mat.gz \
            --plotTitle "H3K27me3 in ${kmean} regions" \
            --samplesLabel "2M" "3M" "10M" "15M" "19M" "24M" \
            --colors "#FF0000" "#FF5500" "#FFAA00" "#00AAFF" "#0055FF" "#0000FF"\
            --plotHeight 10 \
            --plotWidth 12 \
            --regionsLabel "Regions" \
            --yAxisLabel "Signal" \
            --legendLocation "upper-right" \
            --refPointLabel "Center" \
            --perGroup \
            -out  ${data_path}plot/${kmean}_H3K27me3.pdf
    done