antibody=H3K27me3
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)

for tissue in ${tissues[@]}
do
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
    bed=${data_path}all/${antibody}/bed/${antibody}_edd_domain_merged.bed
    saf=${data_path}all/${antibody}/bed/${antibody}_edd_domain_merged.saf
    files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/bed_to_saf.sh ${bed} ${saf}
    featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${antibody}_edd_domain_merged.counts ${files} -F SAF -T 8 
done