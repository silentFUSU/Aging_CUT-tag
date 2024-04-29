tissue=stomach
#创建相关目录，并把bam和bigwig文件转移到目标目录下
bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/mkdir4samples.sh ${tissue}
CUT_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
ATAC_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240402-2_LLX/
# antibodys=(ATAC ATAC H3K9me3 H3K36me3 H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K9me3 H3K36me3 H3K27ac H3K27me3 H3K4me1 H3K4me3 ATAC ATAC H3K9me3 H3K36me3 H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K9me3 H3K36me3 H3K27ac H3K27me3 H3K4me1 H3K4me3)
# antibodys=(ATAC ATAC H3K27ac H3K9me3 H3K27me3 H3K36me3 H3K4me3 H3K4me1 H4K16ac H3K27ac H3K9me3 H3K27me3 H3K36me3 H3K4me3 H3K4me1 H4K16ac)
# samples=$(find ${data_path}bigWig -type f -name "*.bw" -exec basename {} \; | sed 's/_.*//'  | sort) 
antibodys=(ATAC ATAC H3K9me3 H3K36me3 H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K9me3 H3K36me3 H3K27ac H3K27me3 H3K4me1 H3K4me3)
samples=$(find ${data_path}bigWig -type f -name "*.bw" -exec basename {} \; | sed 's/\(^[^_]*_[^_]*\).*$/\1/'  | sort)
# samples=(NTY137_S81_L008 NTY138_S82_L008 NTY139_S83_L008 NTY140_S84_L008 NTY141_S85_L008 NTY142_S86_L008 NTY143_S87_L008 NTY144_S88_L008 NTY145_S89_L008 NTY146_S90_L008 NTY147_S91_L008 NTY148_S92_L008 NTY149_S93_L008 NTY150_S94_L008 NTY151_S95_L008 NTY152_S96_L008)
# samples=$(find ${data_path}bigWig -type f -name "*.bw" -exec basename {} \; | sed 's/\(^[^.]*\).*$/\1/'  | sort)
i=0
for sample in ${samples[@]} 
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

data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240329_LLX_sample/
samples=$(find ${data_path}bigWig -type f -name "*.bw" -exec basename {} \; | sed 's/\(^[^_]*_[^_]*\).*$/\1/'  | sort) 
i=0
antibodys=(ATAC ATAC H3K9me3 H3K36me3 H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K9me3 H3K36me3 H3K27ac H3K27me3 H3K4me1 H3K4me3) 
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
bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/call_peak/bins_calling.sh ${tissue} 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/logs/${tissue}_bin_calling.log &
bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/call_peak/sicer2_mergebam_callpeaks_plot.sh -i /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ \
    -o /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/ -t ${tissue} 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/logs/${tissue}_sicer2_calling.log &
bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/call_peak/macs2_mergebam_callpeaks_plot.sh -i /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ \
    -o /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/ -t ${tissue} 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/logs/${tissue}_macs2_calling.log &

bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/quality_control/FRiP_remove_black_list.sh ${tissue}  2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/logs/${tissue}_FRiP_remove_blacklist.log &