antibodys=(H3K27ac H3K4me1 H3K4me3 ATAC)
# antibodys=(ATAC)
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=$1
ref=mm
for antibody in ${antibodys[@]}
do 
    if [ $antibody = "ATAC" ]; then
        data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
    else
        data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ 
    fi
    bam=($(ls ${data_path}${tissue}/${antibody}/bam/*.bam))
    declare -a young
    declare -a old
    young=()
    old=()
    if [[ $tissue = "brain" && $antibody = "ATAC" ]]; then
        young+=("${bam[0]}")
        young+=("${bam[1]}")
        old+=("${bam[2]}")
        old+=("${bam[3]}")
    elif [ $tissue = "ovary" ]; then
        old+=("${bam[0]}")
        old+=("${bam[2]}")
        young+=("${bam[1]}")
        young+=("${bam[3]}")   
    else
        young+=("${bam[0]}")
        young+=("${bam[2]}")
        old+=("${bam[1]}")
        old+=("${bam[3]}")
    fi
    echo ${young[@]}
    echo ${old[@]}
    echo "samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young[@]} -@ 16"
    echo "samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old[@]} -@ 16"
    
    samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young[@]} -@ 16 &
    samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old[@]} -@ 16 &
    wait
    samtools index ${data_path}${tissue}/${antibody}/tmp.young.merge.bam -@ 16 &
    samtools index ${data_path}${tissue}/${antibody}/tmp.old.merge.bam -@ 16 &
    wait
    mkdir ${data_path}${tissue}/${antibody}/peaks/
    mkdir ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak
    macs2 callpeak -t ${data_path}${tissue}/${antibody}/tmp.young.merge.bam  -f BAMPE -n ${antibody}_young --outdir ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak -g ${ref} --nomodel -q 0.0001  --keep-dup all &
    macs2 callpeak -t ${data_path}${tissue}/${antibody}/tmp.old.merge.bam  -f BAMPE -n ${antibody}_old --outdir ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak -g ${ref} --nomodel -q 0.0001  --keep-dup all &
    wait
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak/${antibody}_young_peaks.narrowPeak |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_narrowpeak.bed  
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak/${antibody}_old_peaks.narrowPeak |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_old_narrowpeak.bed  
    rm ${data_path}${tissue}/${antibody}/tmp*
    bedtools subtract -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_old_narrowpeak.bed -b ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_narrowpeak.bed  > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_old_only_narrowpeak.bed 
    bedtools intersect -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_old_narrowpeak.bed -b ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_narrowpeak.bed  > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_old_intersect_narrowpeak.bed 
    cat ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_narrowpeak.bed  ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_old_narrowpeak.bed | \
        sort -k1,1 -k2,2n | \
        bedtools merge > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_old_narrowpeak.bed
done