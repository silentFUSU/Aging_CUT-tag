# antibodys=(H3K27ac H3K4me1 H3K4me3)
tissue=$1
# for antibody in ${antibodys[@]} 
# do
#     if [ $antibody = ATAC ]; then
#         data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
#     else
#         data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
#     fi
#     samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name *.bam -exec basename {} \; | sed 's/\..*//')
#     echo ${samples[@]}
#     bed=${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_old_narrowpeak.bed
#     saf=${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_old_narrowpeak.saf
#     files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
#     bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${bed} ${saf}
#     featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${antibody}_macs_young_old_narrowpeak.counts ${files} -F SAF -T 8 
# done

antibodys=(H3K27me3 H3K9me3 H3K36me3)
window_size=1000
gap_size=3000
e_value=100
for antibody in ${antibodys[@]}
do
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
    samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name *.bam -exec basename {} \; | sed 's/\..*//')
    echo ${samples[@]}
    bed=${data_path}${tissue}/${antibody}/bed/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed
    saf=${data_path}${tissue}/${antibody}/bed/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}.saf
    files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${bed} ${saf}
    featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}.counts ${files} -F SAF -T 8 
done