# antibodys=(H3K27me3 H3K36me3 H3K9me3)
antibodys=(H3K9me3)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/cellular_aging_GSE133292/Chip_seq/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=$1
ref=hg38
for antibody in ${antibodys[@]}
do 
    search_table=${data_path}CUTTag_search_table_used_in_diff_batch.csv
    cleaned_file=$(mktemp)  
    cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
   #  samples=$(awk -F',' -v t="$tissue" -v y="young" -v a="$antibody" -v c="control" 'NR > 1 && ($1 == t) && ($5 == y) && ($2 == a) && ($4 == c) {print $3}' "$cleaned_file")
   samples=$(awk -F',' -v t="$tissue" -v y="young" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == y) && ($2 == a) {print $3}' "$cleaned_file")
    IFS=$'\n' read -rd '' -a young_array <<<"$samples"  
    echo ${young_array[@]}
    young=()
    for sample in ${young_array[@]}
    do
       file=$(ls ${data_path}bam/${sample}*.bam)
       young+=("$file")
    done
    echo ${young[@]}
   #  samples=$(awk -F',' -v t="$tissue" -v y="old"  -v a="$antibody" -v c="control" 'NR > 1 && ($1 == t) && ($5 == y) && ($2 == a) && ($4 == c) {print $3}' "$cleaned_file")
   samples=$(awk -F',' -v t="$tissue" -v y="old"  -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == y) && ($2 == a) {print $3}' "$cleaned_file")
    IFS=$'\n' read -rd '' -a old_array <<<"$samples"
    echo ${old_array[@]}
    old=()
    for sample in ${old_array[@]}
    do
       file=$(ls ${data_path}bam/${sample}*.bam)
       old+=("$file")
    done
    echo ${old[@]}
    echo "samtools merge -o ${data_path}tmp.young.merge.bam ${young[@]} -@ 16"
    echo "samtools merge -o ${data_path}tmp.old.merge.bam ${old[@]} -@ 16"
    samtools merge -f -o ${data_path}tmp.young.merge.bam ${young[@]} -@ 16 &
    samtools merge -f -o ${data_path}tmp.old.merge.bam ${old[@]} -@ 16 &
    wait
    samtools index ${data_path}tmp.young.merge.bam -@ 16 &
    samtools index ${data_path}tmp.old.merge.bam -@ 16 &
    wait
   #  window_size=1000
   #  gap_size=3000
    window_size=5000
    gap_size=10000
    e_value=100
    mkdir -p ${data_path}peaks/
    sicer  -t ${data_path}tmp.young.merge.bam  -o ${data_path}peaks  -s ${ref} -w ${window_size} -rt 16 -f 300 -egf 0.8 -fdr 0.01 -g ${gap_size} -e ${e_value} -cpu 21 &
    sicer  -t ${data_path}tmp.old.merge.bam  -o ${data_path}peaks  -s ${ref} -w ${window_size} -rt 16 -f 300 -egf 0.8 -fdr 0.01 -g ${gap_size} -e ${e_value} -cpu 21 &
    wait
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR, $4}' ${data_path}peaks/tmp.young.merge-W${window_size}-G${gap_size}.scoreisland |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}bed/${antibody}_young_merge-W${window_size}-G${gap_size}-E${e_value}.bed  
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR, $4}' ${data_path}peaks/tmp.old.merge-W${window_size}-G${gap_size}.scoreisland |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}bed/${antibody}_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed  
    rm ${data_path}tmp*
    bedtools subtract -a ${data_path}bed/${antibody}_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed -b ${data_path}bed/${antibody}_young_merge-W${window_size}-G${gap_size}-E${e_value}.bed > ${data_path}bed/${antibody}_old_only_merge-W${window_size}-G${gap_size}-E${e_value}.bed
    bedtools intersect -a ${data_path}bed/${antibody}_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed -b ${data_path}bed/${antibody}_young_merge-W${window_size}-G${gap_size}-E${e_value}.bed >  ${data_path}bed/${antibody}_young_old_intersect-W${window_size}-G${gap_size}-E${e_value}.bed
    cat  ${data_path}bed/${antibody}_young_merge-W${window_size}-G${gap_size}-E${e_value}.bed ${data_path}bed/${antibody}_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed | \
        sort -k1,1 -k2,2n | \
        bedtools merge  > ${data_path}bed/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed
done