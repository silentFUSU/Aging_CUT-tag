H3K9me3_peaks=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_all_regions.saf
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
for tissue in ${tissues[@]}
do
    files=$(ls ${data_path}${tissue}/bam/*.sorted.bam)
    featureCounts -p -a ${H3K9me3_peaks} -o ${data_path}${tissue}/counts/${tissue}_H3K9me3_peaks_all_regions.counts ${files} -F SAF -T 8 
done

H3K9me3_peaks=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_recursion.saf
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
for tissue in ${tissues[@]}
do
    files=$(ls ${data_path}${tissue}/bam/*.sorted.bam)
    featureCounts -p -a ${H3K9me3_peaks} -o ${data_path}${tissue}/counts/${tissue}_H3K9me3_peaks.counts ${files} -F SAF -T 8 
done