data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=$1
resolution=25000_optimal_parameter
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
cleaned_file=$(mktemp)  
cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
samples=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$cleaned_file")
IFS=$'\n' read -rd '' -a samples_array <<<"$samples"  
for sample in ${samples[@]}
do
   awk 'BEGIN {FS=","; OFS="\t"} NR > 1 {for (i=1; i<=6; i++) gsub(/"/, "", $i); print $1, $2, $3, $4, $5, $6}' ~/projects/Aging_CUT_Tag/data/samples/HiC/${tissue}/loop/HiCCUPS/${sample}_${resolution}/${sample}_${resolution}_loop.csv > ~/projects/Aging_CUT_Tag/data/samples/HiC/${tissue}/loop/HiCCUPS/${sample}_${resolution}/${sample}_${resolution}_loop.bed  
done

young_samples=$(awk -F',' -v t="$tissue" -v y="3M" 'NR > 1 && ($1 == t) && ($5 == y) {print $3}' "$cleaned_file")
old_samples=$(awk -F',' -v t="$tissue" -v y="24M" 'NR > 1 && ($1 == t) && ($5 == y) {print $3}' "$cleaned_file")  
IFS=$'\n' read -rd '' -a young_array <<<"$young_samples"  
IFS=$'\n' read -rd '' -a old_array <<<"$old_samples"  
loop_path=~/projects/Aging_CUT_Tag/data/samples/HiC/${tissue}/loop/HiCCUPS/
merge2Dbed.pl -res ${resolution} ${loop_path}${young_array[0]}_${resolution}/${young_array[0]}_${resolution}_loop.bed ${loop_path}${young_array[1]}_${resolution}/${young_array[1]}_${resolution}_loop.bed \
    ${loop_path}${old_array[0]}_${resolution}/${old_array[0]}_${resolution}_loop.bed ${loop_path}${old_array[1]}_${resolution}/${old_array[1]}_${resolution}_loop.bed -loop > ${loop_path}merged_${resolution}_loop.bed