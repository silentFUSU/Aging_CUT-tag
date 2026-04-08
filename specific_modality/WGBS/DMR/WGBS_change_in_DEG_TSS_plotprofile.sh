tissue=$1
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/${tissue}/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/WGBS/${tissue}/
RNA_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/WGBS_search_table.csv
mkdir -p ${result_path}WGBS_change_in_DEG
mkdir -p ${result_path}WGBS_change_in_DEG/matrix
mkdir -p ${result_path}WGBS_change_in_DEG/plot
cleaned_file=$(mktemp)  
cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
young_samples=$(awk -F',' -v t="$tissue"  'NR > 1 && ($1 == t) && ($5 == "3M") {print $3}' "$cleaned_file") 
old_samples=$(awk -F',' -v t="$tissue"  'NR > 1 && ($1 == t) && ($5 == "24M") {print $3}' "$cleaned_file") 
IFS=$'\n' read -r -d '' -a young_array < <(echo "$young_samples" && printf '\0') 
IFS=$'\n' read -r -d '' -a old_array < <(echo "$old_samples" && printf '\0') 
young1_bw=${data_path}/bw/${young_array[0]}*CpG.bw
young2_bw=${data_path}/bw/${young_array[1]}*CpG.bw
old1_bw=${data_path}/bw/${old_array[0]}*CpG.bw
old2_bw=${data_path}/bw/${old_array[1]}*CpG.bw
/usr/local/lib64/R/bin/Rscript ~/projects/Aging_CUT_Tag/code/WGBS/DMR/WGBS_DEG_TSS_make.R ${tissue}
conditions=(increase decrease)
for condition in ${conditions[@]}
do
    bed=${RNA_path}${tissue}/bed/${condition}_gene_TSS.bed
    computeMatrix reference-point -S $young1_bw $young2_bw $old1_bw $old2_bw -R $bed \
        -a 10000 -b 10000 \
        --numberOfProcessors 10 \
        --skipZeros -o ${result_path}WGBS_change_in_DEG/matrix/WGBS_change_in_DEG_TSS_${condition}.mat.gz

    plotProfile -m ${result_path}WGBS_change_in_DEG/matrix/WGBS_change_in_DEG_TSS_${condition}.mat.gz \
        --plotTitle "WGBS Change in Gene expression ${condition}" \
        --samplesLabel "Young1" "Young2" "Old1" "Old2" \
        --colors "blue" "blue" "red" "red" \
        --plotHeight 10 \
        --plotWidth 12 \
        --regionsLabel "Regions" \
        --yAxisLabel "Signal" \
        --legendLocation "upper-right" \
        --refPointLabel "Center" \
        --perGroup \
        -out ${result_path}WGBS_change_in_DEG/plot/WGBS_change_in_DEG_TSS_${condition}.pdf

    /usr/local/lib64/R/bin/Rscript ~/projects/Aging_CUT_Tag/code/WGBS/DMR/WGBS_change_in_DEG_TSS_plotprofile.R ${tissue} ${condition} &
    
done