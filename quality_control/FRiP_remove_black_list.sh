data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/QC/FRiP/
# tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum bonemarrow ileum heart thymus stomach skin aorta tongue bladder CB jejunum uterus ovary)
# tissues=(BAT)
tissues=(lung)
ref=mm
blacklist=~/ref_data/mm10-blacklist.v2.bed

# antibodys=(H3K4me1)
antibodys=(H3K4me1 H3K4me3 H3K27ac)
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        mkdir ${data_path}${tissue}/${antibody}/FRiP/
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
        mkdir ${data_path}${tissue}/${antibody}/FRiP/bed/

        bedtools intersect -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_old_narrowpeak.bed \
            -b ${blacklist} -v > ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_narrowpeak_rm_blacklist.bed
        bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh \
            ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_narrowpeak_rm_blacklist.bed \
            ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_narrowpeak_rm_blacklist.saf
        
        bedtools intersect -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_old_intersect_narrowpeak.bed \
            -b ${blacklist} -v > ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_intersect_narrowpeak_rm_blacklist.bed
        bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh \
            ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_intersect_narrowpeak_rm_blacklist.bed \
            ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_intersect_narrowpeak_rm_blacklist.saf


        featureCounts -p -a ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_narrowpeak_rm_blacklist.saf \
             -o ${result_path}/union/${tissue}_${antibody}_macs_young_old_narrowpeak_rm_blacklist.counts ${files} -F SAF -T 16 &

        featureCounts -p -a ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_intersect_narrowpeak_rm_blacklist.saf \
             -o ${result_path}/intersect/${tissue}_${antibody}_macs_young_old_intersect_narrowpeak_rm_blacklist.counts ${files} -F SAF -T 16 &
        
        wait
        # rm ${result_path}/union/${tissue}_${antibody}_macs_young_old_narrowpeak_rm_blacklist.counts
        # rm ${result_path}/intersect/${tissue}_${antibody}_macs_young_old_intersect_narrowpeak_rm_blacklist.counts

    done
done

antibodys=(H3K27me3 H3K9me3 H3K36me3) 
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        mkdir ${data_path}${tissue}/${antibody}/FRiP/
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
        mkdir ${data_path}${tissue}/${antibody}/FRiP/bed/
        bedtools intersect -a ${data_path}${tissue}/${antibody}/bed/${antibody}_young_old_merge-W1000-G3000-E100.bed \
            -b ${blacklist} -v > ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_young_old_merge-W1000-G3000-E100_rm_blacklist.bed

        bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh \
            ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_young_old_merge-W1000-G3000-E100_rm_blacklist.bed \
            ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_young_old_merge-W1000-G3000-E100_rm_blacklist.saf
        
        bedtools intersect -a ${data_path}${tissue}/${antibody}/bed/${antibody}_young_old_intersect-W1000-G3000-E100.bed \
            -b ${blacklist} -v > ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_young_old_intersect-W1000-G3000-E100_rm_blacklist.bed

        bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh \
            ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_young_old_intersect-W1000-G3000-E100_rm_blacklist.bed \
            ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_young_old_intersect-W1000-G3000-E100_rm_blacklist.saf


        featureCounts -p -a ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_young_old_merge-W1000-G3000-E100_rm_blacklist.saf \
             -o ${result_path}/union/${tissue}_${antibody}_young_old_merge-W1000-G3000-E100_rm_blacklist.counts ${files} -F SAF -T 16 &

        featureCounts -p -a ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_young_old_intersect-W1000-G3000-E100_rm_blacklist.saf \
             -o ${result_path}/intersect/${tissue}_${antibody}_young_old_intersect-W1000-G3000-E100_rm_blacklist.counts ${files} -F SAF -T 16 &
        
        wait
        # rm ${result_path}/union/${tissue}_${antibody}_young_old_merge-W1000-G3000-E100_rm_blacklist.counts
        # rm ${result_path}/intersect/${tissue}_${antibody}_young_old_intersect-W1000-G3000-E100_rm_blacklist.counts
    done
done

antibodys=(ATAC)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        mkdir ${data_path}${tissue}/${antibody}/FRiP/
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
        mkdir ${data_path}${tissue}/${antibody}/FRiP/bed/

        bedtools intersect -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_old_narrowpeak.bed \
            -b ${blacklist} -v > ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_narrowpeak_rm_blacklist.bed
        bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh \
            ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_narrowpeak_rm_blacklist.bed \
            ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_narrowpeak_rm_blacklist.saf
        
        bedtools intersect -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_old_intersect_narrowpeak.bed \
            -b ${blacklist} -v > ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_intersect_narrowpeak_rm_blacklist.bed
        bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh \
            ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_intersect_narrowpeak_rm_blacklist.bed \
            ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_intersect_narrowpeak_rm_blacklist.saf


        featureCounts -p -a ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_narrowpeak_rm_blacklist.saf \
             -o ${result_path}/union/${tissue}_${antibody}_macs_young_old_narrowpeak_rm_blacklist.counts ${files} -F SAF -T 16 &

        featureCounts -p -a ${data_path}${tissue}/${antibody}/FRiP/bed/${antibody}_macs_young_old_intersect_narrowpeak_rm_blacklist.saf \
             -o ${result_path}/intersect/${tissue}_${antibody}_macs_young_old_intersect_narrowpeak_rm_blacklist.counts ${files} -F SAF -T 16 &
        
        wait
        # rm ${result_path}/union/${tissue}_${antibody}_macs_young_old_narrowpeak_rm_blacklist.counts
        # rm ${result_path}/intersect/${tissue}_${antibody}_macs_young_old_intersect_narrowpeak_rm_blacklist.counts

    done
done