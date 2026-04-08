data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
ref=mm10  
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
antibodys=(H3K27me3 H3K9me3 H3K36me3 H3K27ac H3K4me1 H3K4me3)
for antibody in ${antibodys[@]}
do 
    for tissue in ${tissues[@]}
    do
        # bed=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set/${tissue}.bed
        saf=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/${tissue}/ATAC/bed/ATAC_macs_young_old_narrowpeak_summits_spm3.saf
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
        # bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/bed_to_saf.sh ${bed} ${saf}
        featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${antibody}_in_ATAC_macs_young_old_narrowpeak_summits_spm3.counts ${files[@]} -F SAF -T 8 
    done
done