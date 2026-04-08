H3K9me3_peaks=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_in_or_outpeaks_all_regions.saf
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
# antibodys=(H3K27ac H3K4me3 H3K4me1)
antibodys=(H3K27me3 H3K36me3)
for antibody in ${antibodys[@]}
do
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
    for tissue in ${tissues[@]}
    do
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
        featureCounts -p -a ${H3K9me3_peaks} -o ${data_path}${tissue}/${antibody}/${tissue}_H3K9me3_peaks_all_regions.counts ${files} -F SAF -T 8 
    done
done

H3K9me3_peaks=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.saf
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
for antibody in ${antibodys[@]}
do
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
    for tissue in ${tissues[@]}
    do
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
        featureCounts -p -a ${H3K9me3_peaks} -o ${data_path}${tissue}/${antibody}/${tissue}_H3K9me3_peaks.counts ${files} -F SAF -T 8 
    done
done

H3K9me3_peaks=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_in_or_outpeaks_all_regions.saf
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
antibodys=(ATAC)
for antibody in ${antibodys[@]}
do
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
    for tissue in ${tissues[@]}
    do
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
        featureCounts -p -a ${H3K9me3_peaks} -o ${data_path}${tissue}/${antibody}/${tissue}_H3K9me3_peaks_all_regions.counts ${files} -F SAF -T 8 
    done
done

H3K9me3_peaks=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.saf
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
antibodys=(ATAC)
for antibody in ${antibodys[@]}
do
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
    for tissue in ${tissues[@]}
    do
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
        featureCounts -p -a ${H3K9me3_peaks} -o ${data_path}${tissue}/${antibody}/${tissue}_H3K9me3_peaks.counts ${files} -F SAF -T 8 
    done
done