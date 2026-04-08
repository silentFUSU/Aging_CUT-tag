tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
antibody=H3K9me3

for tissue in ${tissues[@]}
do
    files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}WGBS/${tissue}/hmr/PRC2_LMR_top1000.bed ${data_path}WGBS/${tissue}/hmr/PRC2_LMR_top1000.saf
    featureCounts -p -a ${data_path}WGBS/${tissue}/hmr/PRC2_LMR_top1000.saf -o ${data_path}${tissue}/${antibody}/${tissue}_H3K9me3_PRC2_LMR_top1000.counts ${files} -F SAF -T 8 
done