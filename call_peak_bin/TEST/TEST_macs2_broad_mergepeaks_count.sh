data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=liver
antibody=H3K36me3
ref=mm
bam=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
samtools merge -o ${data_path}${tissue}/${antibody}/tmp.merge.bam ${bam} -@ 16
samtools index ${data_path}${tissue}/${antibody}/tmp.merge.bam -@ 16
mkdir ${data_path}${tissue}/${antibody}/peaks/
mkdir ${data_path}${tissue}/${antibody}/peaks/macs_broadpeak
macs2 callpeak -t ${data_path}${tissue}/${antibody}/tmp.merge.bam  -f BAMPE -n ${antibody} --outdir ${data_path}${tissue}/${antibody}/peaks/macs_broadpeak -g ${ref}  --broad --nomodel --nolambda --keep-dup all
awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' ${data_path}${tissue}/${antibody}/peaks/macs_broadpeak/${antibody}_peaks.broadPeak |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_broad_peak.bed  
bedtools merge -i ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_broad_peak.bed -d 5000 > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_broad_peak_5k_merge.bed
awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_broad_peak_5k_merge.bed > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_broad_peak_5k_merge_2.bed 

# samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.nodup.bam" -exec basename {} \; | sed 's/\..*//')
# for sample in ${samples[@]}
# do
#     bedtools coverage -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_broad_peak_5k_merge_2.bed  -b ${data_path}${tissue}/${antibody}/bam/${sample}.nodup.bam > ${data_path}${tissue}/${antibody}/bed/${sample}_macs_broad_peak_5k_merge_coverage_output.txt  
# done

bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_broad_peak_5k_merge_2.bed  ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_broad_peak_5k_merge_2.saf

files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
featureCounts -p -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_broad_peak_5k_merge_2.saf -o ${data_path}${tissue}/${antibody}/${antibody}_macs_broad_peak_5k_merge.counts ${files} -F SAF -T 8 
