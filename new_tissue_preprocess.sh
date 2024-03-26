tissue=bonemarrow
#创建相关目录，并把bam和bigwig文件转移到目标目录下
bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/mkdir4samples.sh ${tissue}
CUT_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
ATAC_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240301_LLX/
antibodys=(ATAC ATAC H3K9me3 H3K36me3 H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K9me3 H3K36me3 H3K27ac H3K27me3 H3K4me1 H3K4me3)
samples=$(find ${data_path}bigWig -type f -name "*.bw" -exec basename {} \; | sed 's/_.*//'  | sort) 
i=0
for sample in $samples 
do  
    if [ "${antibodys[i]}" = "ATAC" ]; then  
        mv ${data_path}bam/${sample}* ${ATAC_path}${tissue}/${antibodys[i]}/bam/
        mv ${data_path}bigWig/${sample}* ${ATAC_path}${tissue}/${antibodys[i]}/bw/
    else  
        mv ${data_path}bam/${sample}* ${CUT_path}${tissue}/${antibodys[i]}/bam/
        mv ${data_path}bigWig/${sample}* ${CUT_path}${tissue}/${antibodys[i]}/bw/
    fi
    ((i++))
done  

data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240307_LLX/
samples=$(find ${data_path}bigWig -type f -name "*.bw" -exec basename {} \; | sed 's/_.*//'  | sort) 
i=0
  
for sample in $samples 
do  
    if [ "${antibodys[i]}" = "ATAC" ]; then  
        mv ${data_path}bam/${sample}* ${ATAC_path}${tissue}/${antibodys[i]}/bam/
        mv ${data_path}bigWig/${sample}* ${ATAC_path}${tissue}/${antibodys[i]}/bw/
    else  
        mv ${data_path}bam/${sample}* ${CUT_path}${tissue}/${antibodys[i]}/bam/
        mv ${data_path}bigWig/${sample}* ${CUT_path}${tissue}/${antibodys[i]}/bw/
    fi
    ((i++))
done  

#call peak
bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/call_peak/bin_calling.sh ${tissue} 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/logs/${tissue}_bin_calling.log &
bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/call_peak/sicer2_mergebam_callpeaks_plot.sh -i /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ \
    -o /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/ -t ${tissue} 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/logs/${tissue}_sicer2_calling.log &
bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/call_peak/macs2_mergebam_callpeaks_plot.sh -i /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ \
    -o /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/ -t ${tissue} 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/logs/${tissue}_macs2_calling.log &

bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/quality_control/FRiP_remove_black_list.sh ${tissue}  2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/logs/${tissue}_FRiP_remove_blacklist.log &