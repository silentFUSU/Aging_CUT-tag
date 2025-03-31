data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
mkdir -p ${data_path}all/combined_analysis_enhancer/
mkdir -p ${data_path}all/combined_analysis_enhancer/chromHMMbed
antibodys=(H3K27ac H3K27me3 H3K9me3 H3K36me3 H3K4me1 H3K4me3)
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff_batch.csv
ages=(3m 24m)
cleaned_file=$(mktemp)  
cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
> "${data_path}all/combined_analysis_enhancer/cellmarkfiletable.txt"  
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}    
    do
        for age in ${ages[@]}
        do
            if [ "$age" = "3m" ]; then  
                age_label="young"  
            else  
                age_label="old"
            fi  
            samples=$(awk -F',' -v t="$tissue" -v a="$antibody" -v y="$age" 'NR > 1 && ($1 == t) && ($2 == a) && ($5 == y) {print $3}' "$cleaned_file") 
            i=1
            for sample in ${samples[@]}
            do
                echo -e "${tissue}_${age_label}${i}\t${antibody}\t${sample}.bed" >>  ${data_path}all/combined_analysis_enhancer/cellmarkfiletable.txt
                i=$(($i+1))
            done
        done
    done
done