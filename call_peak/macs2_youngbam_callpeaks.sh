antibodys=(H3K27ac H3K4me1 H3K4me3)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ 
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=$1
ref=mm
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

    macs2 callpeak -t ${data_path}${tissue}/${antibody}/tmp.young.merge.bam  -f BAMPE -n ${antibody}_young --outdir ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak -g ${ref} --nomodel -q 0.0001  --keep-dup all &
    macs2 callpeak -t ${data_path}${tissue}/${antibody}/tmp.old.merge.bam  -f BAMPE -n ${antibody}_old --outdir ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak -g ${ref} --nomodel -q 0.0001  --keep-dup all &
    wait
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak/${antibody}_young_peaks.narrowPeak |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_narrowpeak.bed  
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak/${antibody}_old_peaks.narrowPeak |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_old_narrowpeak.bed  
    rm ${data_path}${tissue}/${antibody}/tmp*
    bedtools subtract -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_old_narrowpeak.bed -b ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_narrowpeak.bed  > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_old_only_narrowpeak.bed 
done