data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=$1
resolution=$2
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
cleaned_file=$(mktemp)  
cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
young_samples=$(awk -F',' -v t="$tissue" -v y="3M" 'NR > 1 && ($1 == t) && ($5 == y) {print $3}' "$cleaned_file")
old_samples=$(awk -F',' -v t="$tissue" -v y="24M" 'NR > 1 && ($1 == t) && ($5 == y) {print $3}' "$cleaned_file")  
IFS=$'\n' read -rd '' -a young_array <<<"$young_samples"  
IFS=$'\n' read -rd '' -a old_array <<<"$old_samples"  

TAD_path=${data_path}${tissue}/TAD/homer/
merge2Dbed.pl -res ${resolution} ${TAD_path}${young_array[0]}/${young_array[0]}_${resolution}.tad.2D.bed ${TAD_path}${young_array[1]}/${young_array[1]}_${resolution}.tad.2D.bed \
    ${TAD_path}${old_array[0]}/${old_array[0]}_${resolution}.tad.2D.bed ${TAD_path}${old_array[1]}/${old_array[1]}_${resolution}.tad.2D.bed -tad > ${TAD_path}merged_${resolution}.tad.2D.bed

loop_path=${data_path}${tissue}/loop/homer/
merge2Dbed.pl -res ${resolution} ${loop_path}${young_array[0]}/${young_array[0]}_${resolution}.loop.2D.bed ${loop_path}${young_array[1]}/${young_array[1]}_${resolution}.loop.2D.bed \
    ${loop_path}${old_array[0]}/${old_array[0]}_${resolution}.loop.2D.bed ${loop_path}${old_array[1]}/${old_array[1]}_${resolution}.loop.2D.bed -loop > ${loop_path}merged_${resolution}.loop.2D.bed

result_path=${data_path}${tissue}/differential_analysis
tag_path=${data_path}${tissue}/compartment/homer_compartment/tagDir/
findTADsAndLoops.pl score -tad ${TAD_path}merged_${resolution}.tad.2D.bed -loop ${loop_path}merged_${resolution}.loop.2D.bed \
    -o ${data_path}${tissue}/differential_analysis/homer_score_${resolution} -d ${tag_path}${young_array[0]} ${tag_path}${young_array[1]} ${tag_path}${old_array[0]} ${tag_path}${old_array[1]} \
    -cpu 10 -res ${resolution} 

getDiffExpression.pl  ${data_path}${tissue}/differential_analysis/homer_score_${resolution}.loop.scores.txt young young old old -loop -log2fold 0 > ${data_path}${tissue}/differential_analysis/homer_score_${resolution}.diff.loop.txt
getDiffExpression.pl  ${data_path}${tissue}/differential_analysis/homer_score_${resolution}.tad.scores.txt young young old old -tad -log2fold 0 > ${data_path}${tissue}/differential_analysis/homer_score_${resolution}.diff.tad.txt

findTADsAndLoops.pl score -tad ${TAD_path}merged_${resolution}.tad.2D.bed -loop ${loop_path}merged_${resolution}.loop.2D.bed   \
    -o ${data_path}${tissue}/differential_analysis/homer_score_${resolution}_raw -d ${tag_path}${young_array[0]} ${tag_path}${young_array[1]} ${tag_path}${old_array[0]} ${tag_path}${old_array[1]} \
    -cpu 10 -res ${resolution} 