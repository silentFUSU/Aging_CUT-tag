antibodys=(H3K27me3)
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
for antibody in ${antibodys[@]}
do
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
    bed=${data_path}all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed
    saf=${data_path}all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.saf
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/bed_to_saf.sh ${bed} ${saf}
    for tissue in ${tissues[@]}
    do
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
        featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${tissue}_H3K9me3_peaks_random_200kb_region_rmchrY.counts ${files} -F SAF -T 8 
    done
done

for antibody in ${antibodys[@]}
do
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
    bed=${data_path}all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed
    saf=${data_path}all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.saf
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/bed_to_saf.sh ${bed} ${saf}
    for tissue in ${tissues[@]}
    do
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
        featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${tissue}_H3K9me3_peaks_random_200kb_region_rmchrY_out_recursion_peaks.counts ${files} -F SAF -T 8 
    done
done
