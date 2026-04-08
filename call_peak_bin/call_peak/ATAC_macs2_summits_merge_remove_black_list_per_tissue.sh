tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
for tissue in ${tissues[@]}
do
    # mkdir -p ${data_path}${tissue}/ATAC/peaks/macs_narrowpeak/summits/
    # ln -s ${data_path}${tissue}/ATAC/peaks/macs_narrowpeak/*summits.bed ${data_path}${tissue}/ATAC/peaks/macs_narrowpeak/summits/
    mkdir -p ${data_path}${tissue}/ATAC/peaks/macs_narrowpeak/spm3/
    echo -e "Sample\tGroup\nATAC_young\tyoung\nATAC_old\told" > ${data_path}${tissue}/ATAC/peaks/macs_narrowpeak/file.txt
    Rscript ~/software/ATAC_IterativeOverlapPeakMerging/createIterativeOverlapPeakSet.R \
        --metadata ${data_path}${tissue}/ATAC/peaks/macs_narrowpeak/file.txt \
        --macs2dir ${data_path}${tissue}/ATAC/peaks/macs_narrowpeak/ \
        --outdir ${data_path}${tissue}/ATAC/peaks/macs_narrowpeak/spm3/ \
        --suffix _summits.bed \
        --blacklist ~/ref_data/mm10-blacklist.v2.bed \
        --genome mm10 \
        --spm 3 \
        --rule "(n+1)/2" \
        --extend 250
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}'  ${data_path}${tissue}/ATAC/peaks/macs_narrowpeak/spm3/All_Samples.fwp.filter.non_overlapping.bed |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/ATAC/bed/ATAC_macs_young_old_narrowpeak_summits_spm3.bed
done