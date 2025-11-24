data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=brain
bw=$(ls ${data_path}${tissue}/{H3K27ac,H3K4me,H3K4me3,H3K27me3}/bw/*nodup.bw)
mkdir ${data_path}${tissue}/combined_analysis2_enhancer/matrix
computeMatrix scale-regions -S ${bw} \
  -R ${data_path}${tissue}/combined_analysis2_enhancer/H3K4me1_H3K4me3_logFC15_FDR005_without_promoter.bed    \
  --beforeRegionStartLength 1000 --startLabel start --endLabel end \
  --regionBodyLength 1000 \
  --afterRegionStartLength 1000 \
  --numberOfProcessors 20 \
  --skipZeros -o ${data_path}${tissue}/combined_analysis2_enhancer/matrix/H3K4me1_H3K4me3_logFC15_FDR005_without_promoter.matrix.gz

plotHeatmap -m ${data_path}${tissue}/combined_analysis2_enhancer/matrix/H3K4me1_H3K4me3_logFC15_FDR005_without_promoter.matrix.gz \
     -out ${result_path}${tissue}/deeptools/H3K4me1_H3K4me3_logFC15_FDR005_without_promoter_enhancer.pdf \
     --colorMap Reds \
     --whatToShow 'heatmap and colorbar' --startLabel start --endLabel end \
     --outFileSortedRegions ${data_path}${tissue}/combined_analysis2_enhancer/matrix/H3K4me1_H3K4me3_logFC15_FDR005_without_promoter.bed \
     --kmeans 3