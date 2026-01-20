data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
for tissue in ${tissues[@]}
do
    echo $tissue
    mkdir -p ${data_path}${tissue}/bam/
    search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/WGBS_search_table.csv
    cleaned_file=$(mktemp)  
    cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
    samples=$(awk -F',' -v t="$tissue"  'NR > 1 && ($1 == t) {print $3}' "$cleaned_file")
    IFS=$'\n' read -rd '' -a sample_array <<<"$samples"  
    echo ${sample_array[@]}
    for sample in ${sample_array[@]}
    do
        symlink_path=$(readlink -f ${data_path}${tissue}/bdg/${sample}*_CpG.bdg)
        bed_index=$(echo $symlink_path | grep -b -o '/bed' | cut -d: -f1)
        desired_path=${symlink_path:0:bed_index}
        ln -s ${desired_path}/bam/${sample}_dedup.bam ${data_path}${tissue}/bam/
    done
done