data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
# tissue=$1
# ./your_script.sh "brain liver kidney heart lung spleen"  
ref=mm10  
tissue=(brain liver)
antibodys=(H3K27ac H3K4me3 H3K4me1)
for antibody in ${antibodys[@]}
do 
    bam=()  
    for t in ${tissue[@]}
    do  
        bam+=($(ls ${data_path}${t}/${antibody}/bam/*.bam))  
    done  
    bed=()
    for t in ${tissue[@]}
    do  
        bed+=("$(ls "${data_path}${t}/${antibody}/bed/${antibody}_macs_narrowpeak.bed")")  
    done  
    echo ${bam[@]}
    echo ${bed[@]}
    bedtools intersect -a ${bed[0]} -b ${bed[@]:1}| awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}all/${antibody}/bed/intersect_macs_narrowpeak.bed
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}all/${antibody}/bed/intersect_macs_narrowpeak.bed ${data_path}all/${antibody}/bed/intersect_macs_narrowpeak.saf
    featureCounts -p -a ${data_path}all/${antibody}/bed/intersect_macs_narrowpeak.saf -o ${data_path}all/${antibody}/intersect_macs_narrowpeak.counts ${bam[@]} -F SAF -T 8 

done