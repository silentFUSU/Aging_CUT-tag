tissue=$1
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/${tissue}/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/WGBS/${tissue}/
RNA_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/

mkdir -p ${result_path}WGBS_change_in_DEG
mkdir -p ${result_path}WGBS_change_in_DEG/matrix
mkdir -p ${result_path}WGBS_change_in_DEG/plot
bw=${data_path}bw/${tissue}_Aged_to_Young_1kb_bin_delta_methylation.bw
conditions=(increase decrease)
for condition in ${conditions[@]}
do
    bed=${RNA_path}${tissue}/bed/${condition}_gene_body.bed
    computeMatrix scale-regions -S $bw -R $bed \
        --beforeRegionStartLength 10000 --startLabel TSS --endLabel TES \
        --regionBodyLength 10000 \
        --afterRegionStartLength 10000 \
        --numberOfProcessors 10 \
        --skipZeros -o ${result_path}WGBS_change_in_DEG/matrix/WGBS_1kb_delta_in_DEG_body_${condition}.mat.gz

    plotProfile -m ${result_path}WGBS_change_in_DEG/matrix/WGBS_1kb_delta_in_DEG_body_${condition}.mat.gz \
        --plotTitle "WGBS delta in Gene expression ${condition}" \
        --plotHeight 10 \
        --plotWidth 12 \
        --regionsLabel "Regions" \
        --yAxisLabel "Signal" \
        --legendLocation "upper-right" \
        --refPointLabel "Center" \
        --perGroup \
        -out ${result_path}WGBS_change_in_DEG/plot/WGBS_1kb_delta_in_DEG_body_${condition}.pdf

    /usr/local/lib64/R/bin/Rscript ~/projects/Aging_CUT_Tag/code/WGBS/DMR/WGBS_delta_in_DEG_body_plotprofile.R ${tissue} ${condition} &
done