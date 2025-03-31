data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
ref=mm10
mkdir -p ${data_path}all/combined_analysis_enhancer/
mkdir -p ${data_path}all/combined_analysis_enhancer/chromHMMbed
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff_batch.csv
antibodys=(H3K27ac H3K27me3 H3K9me3 H3K36me3 H3K4me1 H3K4me3)
max_jobs=10
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        samples=$(awk -F',' -v t="$tissue" -v a="$antibody" 'NR > 1 && ($1 == t) && ($2 == a) {print $3}' "$search_table") 
        for sample in ${samples[@]}
        do
            while [ $(current_jobs) -ge $max_jobs ]; do  
                sleep 1  
            done  
            if [ ! -f "${data_path}all/combined_analysis_enhancer/chromHMMbed/${sample}.bed" ]; then  
                echo bedtools bamtobed -i ${data_path}${tissue}/${antibody}/bam/${sample}*.nodup.bam ${data_path}all/combined_analysis_enhancer/chromHMMbed/${sample}.bed
                bedtools bamtobed -i ${data_path}${tissue}/${antibody}/bam/${sample}*.nodup.bam > ${data_path}all/combined_analysis_enhancer/chromHMMbed/${sample}.bed &
            else
                echo ${tissue} ${antibody} ${sample}.bed exist
            fi
        done
    done
done
wait
echo all done