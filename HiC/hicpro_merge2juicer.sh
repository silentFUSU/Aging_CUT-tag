tissues=(brain spleen CB)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
HiC_Pro=/storage/zhangyanxiaoLab/suzhuojie/software/HiC-Pro_3.1.0/
chromsize=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.chrom.sizes
juicer_tools=/storage/zhangyanxiaoLab/suzhuojie/software/juicer/scripts/common/juicer_tools.jar
for tissue in ${tissues[@]}
do
    search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
    cleaned_file=$(mktemp)  
    cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
    samples=$(awk -F',' -v t="$tissue" -v y="3M"  'NR > 1 && ($1 == t) && ($5 == y) {print $3}' "$cleaned_file")
    IFS=$'\n' read -rd '' -a young_array <<<"$samples"  
    echo ${young_array[@]}
    young=()
    for sample in ${young_array[@]}
    do
        file=$(ls ${data_path}${tissue}/ValidPairs/${sample}*.allValidPairs)
        young+=("$file")
    done
    echo ${young[@]}
    samples=$(awk -F',' -v t="$tissue" -v y="24M"  'NR > 1 && ($1 == t) && ($5 == y) {print $3}' "$cleaned_file")
    IFS=$'\n' read -rd '' -a old_array <<<"$samples"
    echo ${old_array[@]}
    old=()
    for sample in ${old_array[@]}
    do
        file=$(ls ${data_path}${tissue}/ValidPairs/${sample}*.allValidPairs)
        old+=("$file")
    done
    echo ${old[@]}

    cat ${young[@]} > ${data_path}${tissue}/ValidPairs/young_combined.allValidPairs &
    cat ${old[@]} > ${data_path}${tissue}/ValidPairs/old_combined.allValidPairs &
    wait
    echo ${tissue}
    bash ${HiC_Pro}bin/utils/hicpro2juicebox.sh -i ${data_path}${tissue}/ValidPairs/young_combined.allValidPairs -g ${chromsize} -j ${juicer_tools} -o ${data_path}${tissue}/juicer/ -t ${data_path}${tissue}/tmp &
    bash ${HiC_Pro}bin/utils/hicpro2juicebox.sh -i ${data_path}${tissue}/ValidPairs/old_combined.allValidPairs -g ${chromsize} -j ${juicer_tools} -o ${data_path}${tissue}/juicer/ -t ${data_path}${tissue}/tmp &
    wait
done