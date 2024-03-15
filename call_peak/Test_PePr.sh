data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=brain
antibody=H3K36me3
young1=LLX146
young2=LLX158
old1=LLX152
old2=LLX164
window_size=default
mkdir ${data_path}${tissue}/${antibody}/PePr/
mkdir ${data_path}${tissue}/${antibody}/PePr/window_size_${window_size}
PePr -c ${data_path}${tissue}/${antibody}/bam/LLX152_S10_L003.nodup.bam,${data_path}${tissue}/${antibody}/bam/LLX164_S106_L003.nodup.bam \
    --chip2 ${data_path}${tissue}/${antibody}/bam/LLX146_S9_L004.nodup.bam,${data_path}${tissue}/${antibody}/bam/LLX158_S103_L003.nodup.bam \
    --output-directory ${data_path}${tissue}/${antibody}/PePr/window_size_${window_size}_diff --windowsize ${window_size} --num-processors 8  \
    --file-format bam --diff

cat ${data_path}${tissue}/${antibody}/PePr/window_size_${window_size}_diff/*chip1_peaks.bed | cut -f 1-3 |sort -k1,1 -k2,2n -|bedtools merge -i stdin | awk -v OFS='\t' 'BEGIN{num=0}{num++; print $0, "peak"num}' - |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/npc_yulab/code/bulk_atac/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/PePr/window_size_${window_size}_diff/old_increase.bed
bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}${tissue}/${antibody}/PePr/window_size_${window_size}_diff/old_increase.bed ${data_path}${tissue}/${antibody}/PePr/window_size_${window_size}_diff/old_increase.saf
files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
featureCounts -p -a ${data_path}${tissue}/${antibody}/PePr/window_size_${window_size}_diff/old_increase.saf -o ${data_path}${tissue}/${antibody}/${antibody}_PePr_old_increase_mergecounts.counts ${files} -F SAF -T 8 

