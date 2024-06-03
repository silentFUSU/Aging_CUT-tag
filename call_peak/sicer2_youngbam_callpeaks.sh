antibodys=(H3K27me3 H3K36me3 H3K9me3)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ 
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=$1
ref=mm10
for antibody in ${antibodys[@]}
do 
    bam=($(ls ${data_path}${tissue}/${antibody}/bam/*.bam))
    declare -a young
    declare -a old
    young=()
    old=()
    young+=("${bam[0]}")
    young+=("${bam[2]}")
    old+=("${bam[1]}")
    old+=("${bam[3]}")

    samtools merge -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young} -@ 16 &
    samtools merge -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old} -@ 16 &
    wait
    samtools index ${data_path}${tissue}/${antibody}/tmp.young.merge.bam -@ 16 &
    samtools index ${data_path}${tissue}/${antibody}/tmp.old.merge.bam -@ 16 &
    wait
    window_size=1000
    gap_size=3000
    # window_size=5000
    # gap_size=10000
    e_value=100
    sicer  -t ${data_path}${tissue}/${antibody}/tmp.young.merge.bam  -o ${data_path}${tissue}/${antibody}/peaks  -s ${ref} -w ${window_size} -rt 16 -f 300 -egf 0.8 -fdr 0.01 -g ${gap_size} -e ${e_value} -cpu 21 &
    sicer  -t ${data_path}${tissue}/${antibody}/tmp.old.merge.bam  -o ${data_path}${tissue}/${antibody}/peaks  -s ${ref} -w ${window_size} -rt 16 -f 300 -egf 0.8 -fdr 0.01 -g ${gap_size} -e ${e_value} -cpu 21 &
    wait
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR, $4}' ${data_path}${tissue}/${antibody}/peaks/tmp.young.merge-W${window_size}-G${gap_size}.scoreisland |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/bed/${antibody}_young_merge-W${window_size}-G${gap_size}-E${e_value}.bed  
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR, $4}' ${data_path}${tissue}/${antibody}/peaks/tmp.old.merge-W${window_size}-G${gap_size}.scoreisland |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/bed/${antibody}_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed  
    rm ${data_path}${tissue}/${antibody}/tmp*
    bedtools subtract -a ${data_path}${tissue}/${antibody}/bed/${antibody}_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed -b ${data_path}${tissue}/${antibody}/bed/${antibody}_young_merge-W${window_size}-G${gap_size}-E${e_value}.bed > ${data_path}${tissue}/${antibody}/bed/${antibody}_old_only_merge-W${window_size}-G${gap_size}-E${e_value}.bed
done