# antibodys=(H3K27ac H3K4me1 H3K4me3 ATAC)
antibodys=(ATAC)
tissue=$1
for antibody in ${antibodys[@]} 
do
    if [ $antibody = "ATAC" ]; then
        data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
    else
        data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
    fi
    samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
    echo ${samples[@]}
    bed=${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_old_narrowpeak.bed
    saf=${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_old_narrowpeak.saf
    files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${bed} ${saf}
    featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${antibody}_macs_young_old_narrowpeak.counts ${files} -F SAF -T 8 
done