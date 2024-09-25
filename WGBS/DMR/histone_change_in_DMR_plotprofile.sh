tissue=liver
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240905_DYQ_005-018_WGBS/DSS_table/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/WGBS/${tissue}/
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table.csv
CUTTag_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
mkdir -p ${result_path}histone_change_in_DMR
mkdir -p ${result_path}histone_change_in_DMR/matrix
mkdir -p ${result_path}histone_change_in_DMR/plot
if [ "$tissue" == "mammarygland" ]; then  
    tissue_label="Mammary Gland"  
else  
    tissue_label=$(echo "$tissue" | awk '{print toupper(substr($0,1,1)) tolower(substr($0,2))}')  
fi 
antibodys=(H3K9me3 H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K27ac)

for antibody in ${antibodys[@]}
do
    young_samples=$(awk -F',' -v t="$tissue_label" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == "3m") && ($2 == a) {print $3}' "$search_table")  
    old_samples=$(awk -F',' -v t="$tissue_label" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == "24m") && ($2 == a) {print $3}' "$search_table")  
    IFS=$'\n' read -r -d '' -a young_array < <(echo "$young_samples" && printf '\0') 
    IFS=$'\n' read -r -d '' -a old_array < <(echo "$old_samples" && printf '\0') 

    # increase_bed=${data_path}bed/DMR_increase.bed
    # decrease_bed=${data_path}bed/DMR_decrease.bed
    young1_bw=${CUTTag_path}${tissue}/${antibody}/bw/${young_array[0]}*.nodup.bw
    young2_bw=${CUTTag_path}${tissue}/${antibody}/bw/${young_array[1]}*.nodup.bw
    old1_bw=${CUTTag_path}${tissue}/${antibody}/bw/${old_array[0]}*.nodup.bw
    old2_bw=${CUTTag_path}${tissue}/${antibody}/bw/${old_array[1]}*.nodup.bw

    conditions=(increase decrease)
    for condition in ${conditions[@]}
    do
        bed=${data_path}bed/${tissue}_DMR_${condition}.bed
        computeMatrix scale-regions -S $young1_bw $old1_bw -R $bed \
            --beforeRegionStartLength 1000 --startLabel start --endLabel end \
            --regionBodyLength 1000 \
            --afterRegionStartLength 1000 \
            --numberOfProcessors 20 \
            --skipZeros -o ${result_path}/histone_change_in_DMR/matrix/${antibody}_rep1_in_WGBS_${condition}.mat.gz 

        computeMatrix scale-regions -S $young2_bw $old2_bw -R $bed \
            --beforeRegionStartLength 1000 --startLabel start --endLabel end \
            --regionBodyLength 1000 \
            --afterRegionStartLength 1000 \
            --numberOfProcessors 20 \
            --skipZeros -o ${result_path}/histone_change_in_DMR/matrix/${antibody}_rep2_in_WGBS_${condition}.mat.gz 

        plotProfile -m ${result_path}/histone_change_in_DMR/matrix/${antibody}_rep1_in_WGBS_${condition}.mat.gz \
            -o ${result_path}/histone_change_in_DMR/plot/${antibody}_rep1_in_WGBS_${condition}.pdf --plotHeight 8 --plotWidth 12 \
            --startLabel start --endLabel end --perGroup \
            --legendLocation best --plotTitle ${tissue}_${antibody}_in_WGBS_${condition}_region

        plotProfile -m ${result_path}/histone_change_in_DMR/matrix/${antibody}_rep2_in_WGBS_${condition}.mat.gz \
            -o ${result_path}/histone_change_in_DMR/plot/${antibody}_rep2_in_WGBS_${condition}.pdf --plotHeight 8 --plotWidth 12 \
            --startLabel start --endLabel end --perGroup \
            --legendLocation best --plotTitle ${tissue}_${antibody}_in_WGBS_${condition}_region
    done
done
/storage/zhangyanxiaoLab/suzhuojie/software/pdfjam-3.11/bin/pdfjam ${result_path}histone_change_in_DMR/plot/{H3K27me3,H3K9me3,H3K36me3,H3K27ac,H3K4me1,H3K4me3}*{increase,decrease}.pdf --nup 4x3  --landscape --outfile ${result_path}histone_change_in_DMR/plot/all_markers.pdf

