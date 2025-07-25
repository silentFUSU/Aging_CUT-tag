bam=()
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/RNA_search_table.csv
cleaned_file=$(mktemp)  
cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
samples=$(awk -F',' 'NR > 1  {print $3}' "$cleaned_file") 
IFS=$'\n' read -r -d '' -a samples_array < <(echo "$samples" && printf '\0') 
for sample in ${samples_array[@]}
do
    file=$(ls ${data_path}*/bam/${sample}*.sorted.bam)
    bam+=("$file")
done
featureCounts -a /storage/zhangyanxiaoLab/share/gtf/mm10.gencode.vM25.annotation.gtf -p -o ${data_path}/all_tissues_combined-chrM.counts ${bam[@]} -F GTF -T 10 -t exon -g gene_name 
