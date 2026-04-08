antibody=H3K27me3
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ 
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=$1
ref=mm10
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff.csv
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

echo "samtools merge -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young[@]} -@ 16"
echo "samtools merge -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old[@]} -@ 16"
samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young[@]} -@ 16 &
samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old[@]} -@ 16 &
wait
samtools index ${data_path}${tissue}/${antibody}/tmp.young.merge.bam -@ 16 &
samtools index ${data_path}${tissue}/${antibody}/tmp.old.merge.bam -@ 16 &
wait
window_size=5000
gap_size=10000
e_value=100
# mkdir -p  ${data_path}${tissue}/${antibody}/peaks/sicer_df/
# sicer_df --treatment_file ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${data_path}${tissue}/${antibody}/tmp.young.merge.bam \
#          --output_directory ${data_path}${tissue}/${antibody}/peaks/sicer_df/ \
#          -s ${ref} -w ${window_size} -rt 16 -f 300 -egf 0.8 -fdr_df 0.05 -fdr 0.05 -g ${gap_size} -e ${e_value} -cpu 21
# awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR, $4}' ${data_path}${tissue}/${antibody}/peaks/sicer_df/tmp.old.merge-W5000-G10000-increased-islands-summary-FDR* |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/peaks/sicer_df/old.merge-W5000-G10000-increased-islands.bed
# awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR, $4}' ${data_path}${tissue}/${antibody}/peaks/sicer_df/tmp.old.merge-W5000-G10000-decreased-islands-summary-FDR* |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/peaks/sicer_df/old.merge-W5000-G10000-decreased-islands.bed
# rm ${data_path}${tissue}/${antibody}/tmp*
mkdir -p  ${data_path}${tissue}/${antibody}/peaks/sicer_control/
sicer  -t ${data_path}${tissue}/${antibody}/tmp.old.merge.bam  -c ${data_path}${tissue}/${antibody}/tmp.young.merge.bam -o ${data_path}${tissue}/${antibody}/peaks/sicer_control/  -s ${ref} -w ${window_size} -rt 16 -f 150 -egf 0.8 -fdr 0.05 -g ${gap_size} -e ${e_value} -cpu 21 &