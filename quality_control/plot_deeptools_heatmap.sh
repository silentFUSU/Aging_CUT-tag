data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
antibodys=(H3K27me3 H3K9me3 H3K36me3 H3K27ac H3K4me3 H3K4me1) 
# antibodys=(H3K27me3 H3K9me3) 
ref=mm10
tissues=(brain liver testis colon kidney)
window_size=1000
gap_size=3000
# window_size=5000
# gap_size=10000
e_value=100
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do 
        bw=($(ls ${data_path}${tissue}/${antibody}/bw/*nodup.bw))  
        temp=${bw[1]}  
        bw[1]=${bw[2]}  
        bw[2]=$temp     
        computeMatrix reference-point -S ${bw[@]} \
            --regionsFileName ${data_path}${tissue}/${antibody}/bed/*_up.bed ${data_path}${tissue}/${antibody}/bed/*_down.bed \
            -a 10000 -b 10000 -bs 20 -p 20 --referencePoint center \
            --skipZeros -o ${data_path}${tissue}/${antibody}/matrix/increase_decrease.mat.gz

        plotHeatmap -m  ${data_path}${tissue}/${antibody}/matrix/increase_decrease.mat.gz \
            -out ${result_path}${tissue}/deeptools/${antibody}_increase_decrease_heatmap.pdf --refPointLabel center \
            --legendLocation best --zMin -5 --zMax 5  --whatToShow 'heatmap and colorbar' \
            --plotTitle "${antibody}"
    done
done
