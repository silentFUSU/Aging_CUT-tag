antibodys=(H3K27me3 H3K9me3 H3K36me3 H3K27ac H3K4me1 H3K4me3 ATAC)
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        if [ $antibody = ATAC ]; then
            data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
        else
            data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
        fi
        samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
        bed=/storage/zhangyanxiaoLab/suzhuojie/ref_data/TE_reference/mm10_TE_all_regions.bed
        saf=/storage/zhangyanxiaoLab/suzhuojie/ref_data/TE_reference/mm10_TE_all_regions.saf
        bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${bed} ${saf}
        featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${antibody}_TE_all.counts ${files} -F SAF -T 8 
    done
done