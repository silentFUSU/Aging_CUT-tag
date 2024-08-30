data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
tissue=ileum
bw=$(ls ${data_path}${tissue}/{H3K27me3,H3K9me3,H3K36me3,H3K27ac,H3K4me1,H3K4me3}/bw/*nodup.bw)
antibodys=(H3K27me3 H3K9me3 H3K36me3)
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
mkdir ${result_path}${tissue}/all_markers_in_others_peak_region/
mkdir ${result_path}${tissue}/all_markers_in_others_peak_region/matrix/
mkdir ${result_path}${tissue}/all_markers_in_others_peak_region/plot/

for antibody in ${antibodys[@]}
do
    mkdir ${result_path}${tissue}/all_markers_in_others_peak_region/matrix/${antibody}
    mkdir ${result_path}${tissue}/all_markers_in_others_peak_region/plot/${antibody}
    bed=${data_path}${tissue}/${antibody}/bed/${antibody}_10kb_in_young_old_merge-W1000-G3000-E100.bed
    computeMatrix scale-regions -S ${bw[@]} \
        -R ${bed} \
        --beforeRegionStartLength 10000 --startLabel start --endLabel end \
        --regionBodyLength 10000 \
        --afterRegionStartLength 10000 \
        --numberOfProcessors 20 \
        --skipZeros -o ${result_path}${tissue}/all_markers_in_others_peak_region/matrix/${antibody}/all_markers_in_${antibody}.mat.gz 
    plotProfile -m  ${result_path}${tissue}/all_markers_in_others_peak_region/matrix/${antibody}/all_markers_in_${antibody}.mat.gz  \
            -out ${result_path}${tissue}/all_markers_in_others_peak_region/plot/${antibody}/all_markers_in_${antibody}.pdf --startLabel start --endLabel end \
            --numPlotsPerRow 4  --legendLocation best \
            --plotTitle "${antibody}"
done
