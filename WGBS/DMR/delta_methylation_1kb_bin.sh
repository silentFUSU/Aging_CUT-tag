data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/
tissue=$1
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/WGBS_search_table.csv
cleaned_file=$(mktemp)  
cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
young_samples_str=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) && ($5 == "3M") {printf "%s_CpG.bdg ", $3}' "$cleaned_file" | xargs) 
sample_files=""  
for sample in $young_samples_str; do  
    sample_files+="${data_path}${tissue}/bdg/$sample "  
done  
cat $sample_files | sort -k1,1 -k2,2n -S 20G > ${data_path}${tissue}/bdg/tmp.young.bdg

old_samples_str=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) && ($5 == "24M") {printf "%s_CpG.bdg ", $3}' "$cleaned_file" | xargs) 
sample_files=""  
for sample in $old_samples_str; do  
    sample_files+="${data_path}${tissue}/bdg/$sample "  
done  
cat $sample_files | sort -k1,1 -k2,2n -S 20G > ${data_path}${tissue}/bdg/tmp.old.bdg

awk '{print $1, $2, $3, $6, $4, $5}' OFS='\t' ${data_path}${tissue}/bdg/tmp.young.bdg | bedtools groupby -i - -g 1,2,3,4 -c 5,6 -o sum | awk '{print $1, $2, $3, $5, $6, $4}' OFS='\t' > ${data_path}${tissue}/bdg/young.bdg
awk '{print $1, $2, $3, $6, $4, $5}' OFS='\t' ${data_path}${tissue}/bdg/tmp.old.bdg | bedtools groupby -i - -g 1,2,3,4 -c 5,6 -o sum | awk '{print $1, $2, $3, $5, $6, $4}' OFS='\t' > ${data_path}${tissue}/bdg/old.bdg

awk '{if($5>4) print $1, $2, $3, 100*$4/$5, $5, $6}' OFS='\t' ${data_path}${tissue}/bdg/young.bdg | bedtools intersect -sorted -wa -wb -a /storage/zhangyanxiaoLab/xiongxiong/Reference/refBed/binGenome/mm10_bin_1k.bed  -b - | awk '{print $1, $2, $3, $7}' OFS='\t' | sort -k1,1 -k2,2n -k3,3n -S 10G > ${data_path}${tissue}/bdg/tmp.young.bed
bedtools groupby -g 1,2,3 -c 4 -o collapse -i ${data_path}${tissue}/bdg/tmp.young.bed | awk '{split($4,a,",")} {if(length(a)>4) print $1, $2, $3}' OFS='\t' > ${data_path}${tissue}/bdg/tmp.young.C5.bin.bed
bedtools intersect -sorted -wa -wb -a ${data_path}${tissue}/bdg/tmp.young.C5.bin.bed -b ${data_path}${tissue}/bdg/tmp.young.bed | cut -f 1-3,7 | bedtools groupby -g 1,2,3 -c 4 -o mean > ${data_path}${tissue}/bdg/young_1kb_bin_CpG.bdg

awk '{if($5>4) print $1, $2, $3, 100*$4/$5, $5, $6}' OFS='\t' ${data_path}${tissue}/bdg/old.bdg | bedtools intersect -sorted -wa -wb -a /storage/zhangyanxiaoLab/xiongxiong/Reference/refBed/binGenome/mm10_bin_1k.bed  -b - | awk '{print $1, $2, $3, $7}' OFS='\t' | sort -k1,1 -k2,2n -k3,3n -S 10G > ${data_path}${tissue}/bdg/tmp.old.bed
bedtools groupby -g 1,2,3 -c 4 -o collapse -i ${data_path}${tissue}/bdg/tmp.old.bed | awk '{split($4,a,",")} {if(length(a)>4) print $1, $2, $3}' OFS='\t' > ${data_path}${tissue}/bdg/tmp.old.C5.bin.bed
bedtools intersect -sorted -wa -wb -a ${data_path}${tissue}/bdg/tmp.old.C5.bin.bed -b ${data_path}${tissue}/bdg/tmp.old.bed | cut -f 1-3,7 | bedtools groupby -g 1,2,3 -c 4 -o mean > ${data_path}${tissue}/bdg/old_1kb_bin_CpG.bdg

bedtools intersect -sorted -wa -wb -a ${data_path}${tissue}/bdg/old_1kb_bin_CpG.bdg -b ${data_path}${tissue}/bdg/young_1kb_bin_CpG.bdg > ${data_path}${tissue}/bdg/tmp.1kb.methylation.bed
awk '{print $1, $2, $3, $4-$8}' OFS='\t' ${data_path}${tissue}/bdg/tmp.1kb.methylation.bed > ${data_path}${tissue}/bdg/${tissue}_Aged_to_Young_1kb_bin_delta_methylation.bdg
bedGraphToBigWig ${data_path}${tissue}/bdg/${tissue}_Aged_to_Young_1kb_bin_delta_methylation.bdg /storage/zhangyanxiaoLab/xiongxiong/index/bismark/mm10/mm10.fa.fai ${data_path}${tissue}/bw/${tissue}_Aged_to_Young_1kb_bin_delta_methylation.bw
rm ${data_path}${tissue}/bdg/tmp*