tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/all/ATAC/macs2_summit_all_tissues_merge/
mkdir -p ${result_path}spm3/
echo -e "Sample\tGroup" > ${result_path}spm3/file.txt
for tissue in ${tissues[@]}
do
    ln -s ${data_path}${tissue}/ATAC/peaks/macs_narrowpeak/ATAC_young_summits.bed ${result_path}spm3/${tissue}_ATAC_young_summits.bed
    ln -s ${data_path}${tissue}/ATAC/peaks/macs_narrowpeak/ATAC_old_summits.bed ${result_path}spm3/${tissue}_ATAC_old_summits.bed
    echo -e "${tissue}_ATAC_young\t${tissue}_young\n${tissue}_ATAC_old\t${tissue}_old" >> ${result_path}spm3/file.txt
done

Rscript ~/software/ATAC_IterativeOverlapPeakMerging/createIterativeOverlapPeakSet.R \
        --metadata ${result_path}spm3/file.txt \
        --macs2dir ${result_path}spm3/ \
        --outdir ${result_path}spm3/ \
        --suffix _summits.bed \
        --blacklist ~/ref_data/mm10-blacklist.v2.bed \
        --genome mm10 \
        --spm 3 \
        --rule "(n+1)/2" \
        --extend 250

awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}'  ${result_path}spm3/All_Samples.fwp.filter.non_overlapping.bed |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${result_path}spm3/ATAC_macs_young_old_narrowpeak_summits_spm3.bed