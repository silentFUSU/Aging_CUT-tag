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
bamCompare -b1 ${data_path}${tissue}/${antibody}/tmp.old.merge.bam \
           -b2 ${data_path}${tissue}/${antibody}/tmp.young.merge.bam \
           --scaleFactorsMethod None --operation log2 --normalizeUsing RPKM \
           -p 10 --binSize 10000 \
           -o ${data_path}${tissue}/${antibody}/bw/old_young_log2ratio.bw