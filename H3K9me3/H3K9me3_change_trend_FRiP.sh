data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
tissue=$1
peaks=${data_path}${tissue}/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_filter.bed
saf=${data_path}${tissue}/H3K9me3/bed/H3K9me3_young_old_merge-W5000-G10000-E100_filter.saf
bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh \
    ${peaks} \
    ${saf}
files=$(ls ${data_path}${tissue}/H3K9me3/bam/*.nodup.bam)
featureCounts -p -a ${saf} \
        -o ${data_path}${tissue}/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_filter.counts ${files} -F SAF -T 16 &
