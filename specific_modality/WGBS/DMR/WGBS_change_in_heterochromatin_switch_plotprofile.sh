tissue=$1
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/${tissue}/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/WGBS/${tissue}/
CUTTag_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/WGBS_search_table.csv
mkdir -p ${result_path}WGBS_change_in_heterochromatin_switch
mkdir -p ${result_path}WGBS_change_in_heterochromatin_switch/matrix
mkdir -p ${result_path}WGBS_change_in_heterochromatin_switch/plot
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
bed=${CUTTag_path}${tissue}/H3K27me3_H3K9me3_intersect/10kb_all_significant_second_quadrant.bed
computeMatrix scale-regions -S $young1_bw $young2_bw $old1_bw $old2_bw -R $bed \
    --beforeRegionStartLength 10000 --startLabel start --endLabel end \
    --regionBodyLength 10000 \
    --afterRegionStartLength 10000 \
    --numberOfProcessors 10 \
    --skipZeros -o ${result_path}WGBS_change_in_heterochromatin_switch/matrix/heterochromatin_switch_WGBS_change.mat.gz

plotProfile -m ${result_path}WGBS_change_in_heterochromatin_switch/matrix/heterochromatin_switch_WGBS_change.mat.gz \
    --plotTitle "WGBS Change in Heterochromatin Switch" \
    --samplesLabel "Young1" "Young2" "Old1" "Old2" \
    --colors "blue" "blue" "red" "red" \
    --plotHeight 10 \
    --plotWidth 12 \
    --regionsLabel "Regions" \
    --yAxisLabel "Signal" \
    --legendLocation "upper-right" \
    --refPointLabel "Center" \
    --perGroup \
    -out ${result_path}WGBS_change_in_heterochromatin_switch/plot/heterochromatin_switch_profile.pdf