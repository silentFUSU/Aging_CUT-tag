tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum bonemarrow ileum heart thymus stomach skin aorta tongue bladder CB jejunum uterus ovary BAT iWAT mammarygland)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/H3K27me3_H3K9me3/H3K27me3_change_in_H3K9me3_peaks/
for tissue in ${tissues[@]}
do
    search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff.csv
    mkdir -p ${result_path}matrix
    mkdir -p ${result_path}plot
    cleaned_file=$(mktemp)  
    cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
    young_samples=$(awk -F',' -v t="$tissue" -v a="H3K27me3" 'NR > 1 && ($1 == t) && ($5 == "3m") && ($2 == a) {print $3}' "$cleaned_file") 
    old_samples=$(awk -F',' -v t="$tissue" -v a="H3K27me3" 'NR > 1 && ($1 == t) && ($5 == "24m") && ($2 == a) {print $3}' "$cleaned_file") 
    IFS=$'\n' read -r -d '' -a young_array < <(echo "$young_samples" && printf '\0') 
    IFS=$'\n' read -r -d '' -a old_array < <(echo "$old_samples" && printf '\0') 
    young1_bw=${data_path}${tissue}/H3K27me3/bw/${young_array[0]}*nodup.bw
    young2_bw=${data_path}${tissue}/H3K27me3/bw/${young_array[1]}*nodup.bw
    old1_bw=${data_path}${tissue}/H3K27me3/bw/${old_array[0]}*nodup.bw
    old2_bw=${data_path}${tissue}/H3K27me3/bw/${old_array[1]}*nodup.bw
    bed=${data_path}${tissue}/H3K9me3/bed/H3K9me3_young_merge-W1000-G3000-E100_compress.bed
    computeMatrix scale-regions -S $young1_bw $young2_bw $old1_bw $old2_bw -R $bed \
        --beforeRegionStartLength 10000 --startLabel Start --endLabel End \
        --regionBodyLength 10000 \
        --afterRegionStartLength 10000 \
        --numberOfProcessors 10 \
        --skipZeros -o ${result_path}matrix/${tissue}_H3K27me3_change_in_H3K9me3_peaks.mat.gz
    
    plotProfile -m ${result_path}matrix/${tissue}_H3K27me3_change_in_H3K9me3_peaks.mat.gz \
        --plotTitle "${tissue} H3K27me3 change in H3K9me3 peaks" \
        --samplesLabel "Young1" "Young2" "Old1" "Old2" \
        --colors "#f38181" "#ff2e63" "#112d4e" "#3f72af" \
        --plotHeight 10 \
        --plotWidth 12 \
        --regionsLabel "Regions" \
        --yAxisLabel "Signal" \
        --legendLocation "upper-right" \
        --refPointLabel "Center" \
        --perGroup \
        --startLabel Start --endLabel End \
        -out ${result_path}plot/${tissue}_H3K27me3_change_in_H3K9me3_peaks.pdf
done

/storage/zhangyanxiaoLab/suzhuojie/software/pdfjam-3.11/bin/pdfjam ${result_path}/plot/{brain,Hip,CB,liver,kidney,testis,colon,cecum,jejunum,ileum,lung,spleen,muscle,pancreas,bonemarrow,heart,thymus,stomach,skin,aorta,tongue,bladder,uterus,ovary,BAT,iWAT,mammarygland}*.pdf --nup 5x6 --landscape --outfile ${result_path}/plot/all_tissues.pdf
/storage/zhangyanxiaoLab/suzhuojie/software/pdfjam-3.11/bin/pdfjam ${result_path}/plot/{brain,Hip,CB,muscle,aorta,tongue,bonemarrow,heart,thymus,testis,colon,spleen,lung,skin,ovary,mammarygland,BAT,iWAT}*.pdf --nup 3x6 --landscape --outfile ${result_path}/plot/tissues_changed.pdf
/storage/zhangyanxiaoLab/suzhuojie/software/pdfjam-3.11/bin/pdfjam ${result_path}/plot/{liver,kidney,cecum,jejunum,ileum,pancreas,stomach,bladder,uterus}*.pdf --nup 3x3 --landscape --outfile ${result_path}/plot/tissues_unchanged.pdf
