data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
# tissue=$1
# ./your_script.sh "brain liver kidney heart lung spleen"  
ref=mm10  
tissue=(brain liver)
antibodys=(H3K27me3 H3K9me3 H3K36me3)
window_size=1000
gap_size=3000
e_value=100
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
        bed+=("$(ls "${data_path}${t}/${antibody}/bed/${antibody}_merge-W${window_size}-G${gap_size}-E${e_value}.bed")")  
    done  
    bedtools intersect -a ${bed[0]} -b ${bed[@]:1}| awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}all/${antibody}/bed/intersect-W${window_size}-G${gap_size}-E${e_value}.bed
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}all/${antibody}/bed/intersect-W${window_size}-G${gap_size}-E${e_value}.bed ${data_path}all/${antibody}/bed/intersect-W${window_size}-G${gap_size}-E${e_value}.saf
    featureCounts -p -a ${data_path}all/${antibody}/bed/intersect-W${window_size}-G${gap_size}-E${e_value}.saf -o ${data_path}all/${antibody}/intersect-W${window_size}-G${gap_size}-E${e_value}.counts ${bam[@]} -F SAF -T 8 

done