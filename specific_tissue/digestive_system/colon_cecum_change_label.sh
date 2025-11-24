tissue=colon
CUT_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/intestine/
ATAC_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/intestine/ATAC/
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/${tissue}/
data_atac_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/${tissue}/

antibodys=(H3K27me3 H3K9me3 H3K36me3 H3K4me3 H3K4me1 H3K27ac ATAC)
tissues_label=(colon cecum colon colon)
for antibody in ${antibodys[@]}
do
    if [ "${antibody}" = "ATAC" ]; then  
        samples=$(find ${data_atac_path}${antibody}/bw -type f -name "*.nodup.bw" -exec basename {} \; | sed 's/\(_[^_]*\).*/\1/' | sort)
    else
        samples=$(find ${data_path}${antibody}/bw -type f -name "*.nodup.bw" -exec basename {} \; | sed 's/\(_[^_]*\).*/\1/' | sort)
    fi
    i=0
    for sample in ${samples[@]}
        do
            if [ "${antibody}" = "ATAC" ]; then  
                cp ${data_atac_path}${antibody}/bam/${sample}* ${ATAC_path}${tissues_label[i]}/${antibody}/bam/ &
                cp ${data_atac_path}${antibody}/bw/${sample}* ${ATAC_path}${tissues_label[i]}/${antibody}/bw/ &
            else  
                cp ${data_path}${antibody}/bam/${sample}* ${CUT_path}${tissues_label[i]}/${antibody}/bam/ & 
                cp ${data_path}${antibody}/bw/${sample}* ${CUT_path}${tissues_label[i]}/${antibody}/bw/ &
            fi
            ((i++))
        done
done

bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/digestive_system/bins_calling.sh ${tissue} 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/${tissue}_corrected_bins_calling.log &
nohup bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/digestive_system/sicer2_age_split_bam_callpeaks.sh ${tissue} 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/${tissue}_corrected_sicer_peak_calling.log &
nohup bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/digestive_system/macs2_age_split_bam_callpeaks.sh ${tissue} 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/${tissue}_corrected_macs2_peak_calling.log &