tissues=(brain CB kidney liver lung bonemarrow colon heart Hip mammarygland stomach thymus skin muscle cecum ileum)
tissues=(pancreas spleen)
for tissue in ${tissues[@]}
do
    bed=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/${tissue}/compartment/homer_compartment/compartment_change_50000.bed
    saf=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/${tissue}/compartment/homer_compartment/compartment_change_50000.saf
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${bed} ${saf}
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
    files=$(ls ${data_path}${tissue}/bam/*.sorted.bam)
    featureCounts -p -a ${saf} -o ${data_path}${tissue}/counts/${tissue}_compartment_50000.counts ${files} -F SAF -T 8 
done

antibodys=(H3K27ac H3K4me3 H3K4me1 H3K27me3 H3K36me3 H3K9me3)
for antibody in ${antibodys[@]}
do
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
    for tissue in ${tissues[@]}
    do
        saf=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/${tissue}/compartment/homer_compartment/compartment_change_50000.saf
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
        featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${tissue}_compartment_50000.counts ${files} -F SAF -T 8 
    done
done

antibodys=(ATAC)
tissues=(brain CB kidney liver lung bonemarrow colon heart Hip mammarygland stomach thymus skin muscle cecum ileum pancreas spleen)
for antibody in ${antibodys[@]}
do
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
    for tissue in ${tissues[@]}
    do
        saf=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/${tissue}/compartment/homer_compartment/compartment_change_50000.saf
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
        featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${tissue}_compartment_50000.counts ${files} -F SAF -T 8 
    done
done