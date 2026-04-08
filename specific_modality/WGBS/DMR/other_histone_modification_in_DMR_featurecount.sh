data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
ref=mm10  
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
antibodys=(H3K27me3 H3K9me3 H3K36me3 H3K27ac H3K4me1 H3K4me3)
conditions=(increase decrease)
for antibody in ${antibodys[@]}
do 
    for tissue in ${tissues[@]}
    do
        for condition in ${conditions[@]}
        do
            bed=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/${tissue}/DSS_table/bed/${tissue}_DMR_${condition}_delta01.bed
            saf=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/${tissue}/DSS_table/bed/${tissue}_DMR_${condition}_delta01.saf
            files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
            bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/bed_to_saf.sh ${bed} ${saf}
            featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${antibody}_DMR_${condition}_delta01.counts ${files[@]} -F SAF -T 8 
        done
    done
done

data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
antibodys=(ATAC)
conditions=(increase decrease)
for antibody in ${antibodys[@]}
do 
    for tissue in ${tissues[@]}
    do
        for condition in ${conditions[@]}
        do
            bed=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/${tissue}/DSS_table/bed/${tissue}_DMR_${condition}_delta01.bed
            saf=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/${tissue}/DSS_table/bed/${tissue}_DMR_${condition}_delta01.saf
            files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
            bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/bed_to_saf.sh ${bed} ${saf}
            featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${antibody}_DMR_${condition}_delta01.counts ${files[@]} -F SAF -T 8 
        done
    done
done