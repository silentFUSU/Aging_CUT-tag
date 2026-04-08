data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/mnt/transposon2/zhangyanxiaoLab/suzhuojie/project/Aging_CUT_Tag/LMJ_Network_rawdata/
antibody=H3K27ac
tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum bonemarrow ileum heart thymus stomach skin aorta tongue bladder CB jejunum uterus ovary BAT iWAT mammarygland)
bed=/mnt/transposon2/zhangyanxiaoLab/suzhuojie/project/Aging_CUT_Tag/LMJ_Network_rawdata/ATAC_union_peak.bed
saf=/mnt/transposon2/zhangyanxiaoLab/suzhuojie/project/Aging_CUT_Tag/LMJ_Network_rawdata/ATAC_union_peak.saf
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${bed} ${saf}
for tissue in ${tissues[@]}
do
    search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff_batch.csv
    cleaned_file=$(mktemp)  
    cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
    samples=$(awk -F',' -v t="$tissue" -v a="$antibody" 'NR > 1 && ($1 == t) && ($2 == a) {print $3}' "$cleaned_file") 
    IFS=$'\n' read -r -d '' -a samples_array < <(echo "$samples" && printf '\0') 
    bams=()
    for sample in ${samples_array[@]}
    do
        for file in $(find "$data_path${tissue}/${antibody}/bam/" -type f -name "${sample}*.bam"); do
            bams+=("$file")
        done
    done
    featureCounts -p -a ${saf} -o ${result_path}/${antibody}/${tissue}_in_ATAC_peaks.counts ${bams[@]} -F SAF -T 8 
done