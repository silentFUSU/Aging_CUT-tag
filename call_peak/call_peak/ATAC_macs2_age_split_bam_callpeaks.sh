# antibodys=(H3K27ac H3K4me1 H3K4me3)
antibodys=(ATAC)
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=$1
ref=mm
for antibody in ${antibodys[@]}
do 
    if [ $antibody = "ATAC" ]; then
        data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
        search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/ATAC_search_table_batch.csv
    else
        data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
        search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff.csv 
    fi
    cleaned_file=$(mktemp)  
    cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
    samples=$(awk -F',' -v t="$tissue" -v y="3m" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == y) && ($2 == a) {print $3}' "$cleaned_file")
    IFS=$'\n' read -rd '' -a young_array <<<"$samples"  
    echo ${young_array[@]}
    young=()
    for sample in ${young_array[@]}
    do
       file=$(ls ${data_path}${tissue}/${antibody}/bam/${sample}*.bam)
       young+=("$file")
    done
    echo ${young[@]}
    samples=$(awk -F',' -v t="$tissue" -v y="24m" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == y) && ($2 == a) {print $3}' "$cleaned_file")
    IFS=$'\n' read -rd '' -a old_array <<<"$samples"
    echo ${old_array[@]}
    old=()
    for sample in ${old_array[@]}
    do
       file=$(ls ${data_path}${tissue}/${antibody}/bam/${sample}*.bam)
       old+=("$file")
    done
    echo ${old[@]}
    echo "samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young[@]} -@ 16"
    echo "samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old[@]} -@ 16"
    
    samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young[@]} -@ 16 &
    samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old[@]} -@ 16 &
    wait
    samtools index ${data_path}${tissue}/${antibody}/tmp.young.merge.bam -@ 16 &
    samtools index ${data_path}${tissue}/${antibody}/tmp.old.merge.bam -@ 16 &
    wait
    mkdir ${data_path}${tissue}/${antibody}/peaks/
    mkdir ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak_01
    macs2 callpeak -t ${data_path}${tissue}/${antibody}/tmp.young.merge.bam  -f BAMPE -n ${antibody}_young --outdir ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak_01 -g ${ref} --nomodel -q 0.0001  --keep-dup all &
    macs2 callpeak -t ${data_path}${tissue}/${antibody}/tmp.old.merge.bam  -f BAMPE -n ${antibody}_old --outdir ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak_01 -g ${ref} --nomodel -q 0.0001  --keep-dup all &
    wait
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak_01/${antibody}_young_peaks.narrowPeak |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_narrowpeak_01.bed  
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak_01/${antibody}_old_peaks.narrowPeak |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_old_narrowpeak_01.bed  
    rm ${data_path}${tissue}/${antibody}/tmp*
    bedtools subtract -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_old_narrowpeak_01.bed -b ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_narrowpeak_01.bed  > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_old_only_narrowpeak_01.bed 
    bedtools intersect -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_old_narrowpeak_01.bed -b ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_narrowpeak_01.bed  > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_old_intersect_narrowpeak_01.bed 
    cat ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_narrowpeak_01.bed  ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_old_narrowpeak_01.bed | \
        sort -k1,1 -k2,2n | \
        bedtools merge > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_old_narrowpeak_01.bed
done