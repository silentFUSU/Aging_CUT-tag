antibody=H3K27me3
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
for tissue in ${tissues[@]}
do
    samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
    files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
    bed=${data_path}/WGBS/${tissue}/hmr/all_samples_hmr.bed
    saf=${data_path}/WGBS/${tissue}/hmr/all_samples_hmr.saf
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/bed_to_saf.sh ${bed} ${saf}
    featureCounts -p -a ${saf} -o ${data_path}${tissue}/${antibody}/${antibody}_all_samples_hmr.counts ${files} -F SAF -T 8 
done