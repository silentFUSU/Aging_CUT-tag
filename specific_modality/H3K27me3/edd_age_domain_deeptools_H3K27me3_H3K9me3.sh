## merge bam plot
# tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum bonemarrow ileum heart thymus stomach skin aorta tongue bladder CB jejunum uterus ovary BAT iWAT mammarygland)
tissues=(brain)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/H3K27me3_domain/
mkdir -p $result_path
max_jobs=4
current_jobs() {  
    jobs -rp | wc -l  
}  

for tissue in ${tissues[@]}
do

    mkdir -p ${result_path}matrix_merge
    mkdir -p ${result_path}plot_merge
    if [[ "$antibody" == "H3K27me3" || "$antibody" == "H3K36me3" || "$antibody" == "H3K9me3" ]]; then
        binsize=1000    
    else
        binsize=50
    fi

    K27_young_bw=${data_path}${tissue}/H3K27me3/bw/young_bs${binsize}.bw
    K27_old_bw=${data_path}${tissue}/H3K27me3/bw/old_bs${binsize}.bw
    K9_young_bw=${data_path}${tissue}/H3K9me3/bw/young_bs${binsize}.bw
    K9_old_bw=${data_path}${tissue}/H3K9me3/bw/old_bs${binsize}.bw

    bed=${data_path}${tissue}/H3K27me3/peaks/edd/edd_peaks_fdr05.bed
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    computeMatrix scale-regions -S $K27_young_bw $K27_old_bw $K9_young_bw $K9_old_bw -R $bed \
        --beforeRegionStartLength 100000 --startLabel Start --endLabel End \
        --regionBodyLength 100000 \
        --afterRegionStartLength 100000 \
        --numberOfProcessors 10 \
        --skipZeros -o ${result_path}matrix_merge/${tissue}_H3K27me3_H3K9me3_change_in_H3K27me3_domain.mat.gz &

done

wait

for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        # plotProfile -m ${result_path}matrix_merge/${tissue}_${antibody}_change_in_H3K27me3_domain.mat.gz \
        #     --plotTitle "${tissue} ${antibody}" \
        #     --samplesLabel "Young"  "Old" \
        #     --colors "#e64b35" "#3c5488" \
        #     --plotHeight 10 \
        #     --plotWidth 10 \
        #     --regionsLabel "Regions" \
        #     --yAxisLabel "Signal" \
        #     --legendLocation "upper-right" \
        #     --refPointLabel "Center" \
        #     --perGroup \
        #     --startLabel Start --endLabel End \
        #     --yMax 0.5 \
        #     -out ${result_path}plot_merge/${tissue}_${antibody}_change_in_H3K27me3_domain.pdf &
        
        plotHeatmap -m ${result_path}matrix_merge/${tissue}_H3K27me3_H3K9me3_change_in_H3K27me3_domain.mat.gz \
            -out ${result_path}plot_merge/${tissue}_H3K27me3_H3K9me3_change_in_H3K27me3_domain_heatmap.pdf \
            --startLabel Start --endLabel End \
            --zMin 0 \
            --zMax 1 \
            --colorList 'white,#3852B4' \
            --samplesLabel "K27 Young" "K27 Old" "K9 Young" "K9 Old" \
            --whatToShow 'heatmap and colorbar'
    done
done