data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/QC/FRiP/
# tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum)
tissues=$1
ref=mm
blacklist=~/ref_data/mm10-blacklist.v2.bed

antibodys=(H3K27ac H3K4me3 H3K4me1)
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
        bedtools intersect -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_narrowpeak.bed \
            -b ${blacklist} -v > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_narrowpeak_rm_blacklist.bed
        bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh \
            ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_narrowpeak_rm_blacklist.bed \
            ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_narrowpeak_rm_blacklist.saf
        featureCounts -p -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_narrowpeak_rm_blacklist.saf \
             -o ${result_path}/${tissue}_${antibody}_rm_blacklist.counts ${files} -F SAF -T 8 
        rm ${result_path}/${tissue}_${antibody}_rm_blacklist.counts
    done
done

antibodys=(H3K27me3 H3K9me3 H3K36me3) 
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
        bedtools intersect -a ${data_path}${tissue}/${antibody}/bed/${antibody}_merge-W1000-G3000-E100.bed \
            -b ${blacklist} -v > ${data_path}${tissue}/${antibody}/bed/${antibody}_merge-W1000-G3000-E100_rm_blacklist.bed
        bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh \
            ${data_path}${tissue}/${antibody}/bed/${antibody}_merge-W1000-G3000-E100_rm_blacklist.bed \
            ${data_path}${tissue}/${antibody}/bed/${antibody}_merge-W1000-G3000-E100_rm_blacklist.saf
        featureCounts -p -a ${data_path}${tissue}/${antibody}/bed/${antibody}_merge-W1000-G3000-E100_rm_blacklist.saf \
             -o ${result_path}/${tissue}_${antibody}_rm_blacklist.counts ${files} -F SAF -T 16
        rm ${result_path}/${tissue}_${antibody}_rm_blacklist.counts
    done
done