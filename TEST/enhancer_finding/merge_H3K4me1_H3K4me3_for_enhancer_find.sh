data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=brain
peaks=$(ls ${data_path}${tissue}/{H3K4me,H3K4me3}/bed/*mergecounts.bed)
cat $peaks | cut -f 1-3 |sort -k1,1 -k2,2n -|bedtools merge -i stdin | awk -v OFS='\t' 'BEGIN{num=0}{num++; print $0, "peak"num}' - |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/npc_yulab/code/bulk_atac/keep_regular_chroms.r > ${data_path}${tissue}/H3K4me1_H3K4me3_merge.bed
bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}${tissue}/H3K4me1_H3K4me3_merge.bed ${data_path}${tissue}/H3K4me1_H3K4me3_merge.saf
files=$(ls ${data_path}${tissue}/{H3K4me,H3K4me3}/bam/*.nodup.bam)
featureCounts -p -a ${data_path}${tissue}/H3K4me1_H3K4me3_merge.saf -o ${data_path}${tissue}/H3K4me1_H3K4me3_merge.counts ${files} -F SAF -T 8 

bedtools intersect -a ${data_path}${tissue}/H3K27ac/bed/H3K27ac_macs2_mergecounts_without_promoter.bed -b ${data_path}${tissue}/H3K4me1_H3K4me3_logFC1_FDR005.bed >  ${data_path}${tissue}/H3K27ac_intersect_H3K4me1_H3K4me3_logFC1_FDR005.bed
bedtools intersect -a ${data_path}${tissue}/H3K27ac/bed/H3K27ac_macs2_mergecounts_without_promoter.bed -b ${data_path}${tissue}/H3K4me1_H3K4me3_logFC2_FDR005.bed >  ${data_path}${tissue}/H3K27ac_intersect_H3K4me1_H3K4me3_logFC2_FDR005.bed
bedtools intersect -a ${data_path}${tissue}/H3K27ac/bed/H3K27ac_macs2_mergecounts_without_promoter.bed -b ${data_path}${tissue}/H3K4me1_H3K4me3_logFC15_FDR005.bed >  ${data_path}${tissue}/H3K27ac_intersect_H3K4me1_H3K4me3_logFC15_FDR005.bed

bedtools intersect -a ${data_path}${tissue}/H3K27ac_intersect_H3K4me1_H3K4me3_logFC15_FDR005.bed -b ${data_path}${tissue}/H3K27me3/bed/H3K27me3_merge-W1000-G3000-E100.bed -v > ${data_path}${tissue}/positive_enhancer.bed