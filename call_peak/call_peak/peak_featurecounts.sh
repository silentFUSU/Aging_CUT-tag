# antibodys=(H3K27ac H3K4me1 H3K4me3)
tissue=$1
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
antibodys=(ATAC)
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]} 
    do
        if [ $antibody = ATAC ]; then
            data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
            bed=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/all/ATAC/macs2_summit_all_tissues_merge/spm3/ATAC_macs_young_old_narrowpeak_summits_spm3.bed
            saf=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/all/ATAC/macs2_summit_all_tissues_merge/spm3/ATAC_macs_young_old_narrowpeak_summits_spm3.saf
        else
            data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
            bed=
            saf=
        fi
        samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name *.bam -exec basename {} \; | sed 's/\..*//')
        echo ${samples[@]}
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
        bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${bed} ${saf}
        featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${antibody}_macs_young_old_narrowpeak_summits_spm3_all_tissues_merge.counts ${files} -F SAF -T 8 
    done
done

# antibodys=(H3K27me3 H3K9me3 H3K36me3)
antibodys=(H3K9me3)
window_size=5000
gap_size=10000
e_value=100
for antibody in ${antibodys[@]}
do
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
    samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name *.bam -exec basename {} \; | sed 's/\..*//')
    echo ${samples[@]}
    # bed=${data_path}${tissue}/${antibody}/bed/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed
    bed=${data_path}all/${antibody}/bed/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}_recursion.bed
    # saf=${data_path}${tissue}/${antibody}/bed/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}.saf
    saf=${data_path}all/${antibody}/bed/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}_recursion.saf
    files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${bed} ${saf}
    featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}_recursion.counts ${files} -F SAF -T 8 
done