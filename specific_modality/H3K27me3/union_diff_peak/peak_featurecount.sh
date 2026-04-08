antibody=H3K27me3
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
tissue_num=20
for tissue in ${tissues[@]}
do
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
    bed=${data_path}all/${antibody}/bed/${antibody}_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_${tissue_num}_tissues.bed
    saf=${data_path}all/${antibody}/bed/${antibody}_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_${tissue_num}_tissues.saf
    files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/bed_to_saf.sh ${bed} ${saf}
    featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${antibody}_young_merge-W5000-G10000-E100_peak_merged_exist_in_larger_${tissue_num}_tissues.counts ${files} -F SAF -T 8 
done