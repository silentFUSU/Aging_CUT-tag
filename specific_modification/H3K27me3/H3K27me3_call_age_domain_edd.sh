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

# echo "samtools merge -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young[@]} -@ 16"
# echo "samtools merge -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old[@]} -@ 16"
# samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young[@]} -@ 16 &
# samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old[@]} -@ 16 &
# wait
# samtools index ${data_path}${tissue}/${antibody}/tmp.young.merge.bam -@ 16 &
# samtools index ${data_path}${tissue}/${antibody}/tmp.old.merge.bam -@ 16 &
# wait

mkdir -p  ${data_path}${tissue}/${antibody}/peaks/edd/
source /storage/zhangyanxiaoLab/suzhuojie/miniconda3/etc/profile.d/conda.sh
conda activate edd
if [[ "$tissue" == "mammarygland" || "$tissue" == "ovary" || "$tissue" == "uterus" ]]; then
    genome_size="/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/edd_mm10_without_chrY.chrom.sizes"
else
    genome_size="/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/edd_mm10.chrom.sizes"
fi

unalignable_regions=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/edd_unalignable_regions.bed
# edd  -p 16 --fdr 0.05  -n 50000 --gap-penalty 80 \
#     ${genome_size} ${unalignable_regions} ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${data_path}${tissue}/${antibody}/peaks/edd/
    
# mv ${data_path}${tissue}/${antibody}/peaks/edd/edd_peaks.bed ${data_path}${tissue}/${antibody}/peaks/edd/edd_peaks_fdr05.bed
# mv ${data_path}${tissue}/${antibody}/peaks/edd/log.txt ${data_path}${tissue}/${antibody}/peaks/edd/log_fdr05.txt

# edd  -p 16 --fdr 2  -n 50000 --gap-penalty 80 \
#     ${genome_size} ${unalignable_regions} ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${data_path}${tissue}/${antibody}/peaks/edd/
# mv ${data_path}${tissue}/${antibody}/peaks/edd/edd_peaks.bed ${data_path}${tissue}/${antibody}/peaks/edd/edd_peaks_all.bed
# mv ${data_path}${tissue}/${antibody}/peaks/edd/log.txt ${data_path}${tissue}/${antibody}/peaks/edd/log_all.txt

edd  -p 16 --fdr 0.05  -n 50000 --gap-penalty 80 \
    ${genome_size} ${unalignable_regions} ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${data_path}${tissue}/${antibody}/tmp.old.merge.bam  ${data_path}${tissue}/${antibody}/peaks/edd/
    
mv ${data_path}${tissue}/${antibody}/peaks/edd/edd_peaks.bed ${data_path}${tissue}/${antibody}/peaks/edd/young_edd_peaks_fdr05.bed
mv ${data_path}${tissue}/${antibody}/peaks/edd/log.txt ${data_path}${tissue}/${antibody}/peaks/edd/young_log_fdr05.txt

edd  -p 16 --fdr 2  -n 50000 --gap-penalty 80 \
    ${genome_size} ${unalignable_regions} ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${data_path}${tissue}/${antibody}/tmp.old.merge.bam  ${data_path}${tissue}/${antibody}/peaks/edd/
mv ${data_path}${tissue}/${antibody}/peaks/edd/edd_peaks.bed ${data_path}${tissue}/${antibody}/peaks/edd/young_edd_peaks_all.bed
mv ${data_path}${tissue}/${antibody}/peaks/edd/log.txt ${data_path}${tissue}/${antibody}/peaks/edd/young_log_all.txt