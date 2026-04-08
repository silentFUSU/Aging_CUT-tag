data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data_transposon2/20250908_MWY_CUTTAG/
bed=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/lung/H3K27me3/peaks/edd/edd_peaks_fdr05.bed
M2_bw=${data_path}bigWig/MWY922-1*nodup.bw
M3_bw=${data_path}bigWig/MWY922-2*nodup.bw
M10_bw=${data_path}bigWig/MWY922-3*nodup.bw
M15_bw=${data_path}bigWig/MWY922-4*nodup.bw
M19_bw=${data_path}bigWig/MWY922-5*nodup.bw
M24_bw=${data_path}bigWig/MWY922-6*nodup.bw
computeMatrix scale-regions -S $M2_bw $M3_bw $M10_bw $M15_bw $M19_bw $M24_bw -R $bed \
    --beforeRegionStartLength 1000 --startLabel Start --endLabel End \
    --regionBodyLength 1000 \
    --afterRegionStartLength 1000 \
    --numberOfProcessors 10 \
    --skipZeros -o ${data_path}matrix/H3K9me3_in_H3K27me3_domains.mat.gz
plotProfile -m ${data_path}matrix/H3K9me3_in_H3K27me3_domains.mat.gz \
    --plotTitle "H3K9me3 in H3K27me3 domain regions" \
    --samplesLabel "2M" "3M" "10M" "15M" "19M" "24M" \
    --colors "#FF0000" "#FF5500" "#FFAA00" "#00AAFF" "#0055FF" "#0000FF"\
    --plotHeight 10 \
    --plotWidth 12 \
    --regionsLabel "Regions" \
    --yAxisLabel "Signal" \
    --legendLocation "upper-right" \
    --refPointLabel "Center" \
    --perGroup \
    -out  ${data_path}plot/H3K9me3_in_H3K27me3_domains.pdf
