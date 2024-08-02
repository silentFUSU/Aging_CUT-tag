data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=testis
state=E13
peaks=$(ls ${result_path}${tissue}/ChromHMM/15/{young1,young2,old1,old2}_*_segments.bed)
mkdir ${data_path}${tissue}/combined_analysis_enhancer/bed/
mkdir ${data_path}${tissue}/combined_analysis_enhancer/bed/15
awk -v state="${state}" 'BEGIN {OFS="\t"} $4 == state {print $1, $2, $3}' ${peaks} |sort -k1,1 -k2,2n > ${data_path}${tissue}/combined_analysis_enhancer/bed/15/tmp_merged_segments_${state}.bed  
bedtools merge -i ${data_path}${tissue}/combined_analysis_enhancer/bed/15/tmp_merged_segments_${state}.bed  -d 10000 > ${data_path}${tissue}/combined_analysis_enhancer/bed/15/tmp2_merged_segments_${state}.bed 
awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' ${data_path}${tissue}/combined_analysis_enhancer/bed/15/tmp2_merged_segments_${state}.bed  > ${data_path}${tissue}/combined_analysis_enhancer/bed/15/merged_segments_${state}.bed 
rm ${data_path}${tissue}/combined_analysis_enhancer/bed/15/tmp*
bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}${tissue}/combined_analysis_enhancer/bed/15/merged_segments_${state}.bed   ${data_path}${tissue}/combined_analysis_enhancer/bed/15/merged_segments_${state}.saf
files=$(ls ${data_path}${tissue}/H3K27me3/bam/*.nodup.bam)
featureCounts -p -a ${data_path}${tissue}/combined_analysis_enhancer/bed/15/merged_segments_${state}.saf -o ${data_path}${tissue}/H3K27me3/H3K27me3_segments_${state}.counts ${files} -F SAF -T 8 
files=$(ls ${data_path}${tissue}/H3K9me3/bam/*.nodup.bam)
featureCounts -p -a ${data_path}${tissue}/combined_analysis_enhancer/bed/15/merged_segments_${state}.saf -o ${data_path}${tissue}/H3K9me3/H3K9me3_segments_${state}.counts ${files} -F SAF -T 8 
