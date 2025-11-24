# conda activate py27
# data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/FA_vs_noFA/
# bash sicer2_mergebam_callpeaks_plot.sh -i /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ -o /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/ -t spleen 
func() {
    echo "Usage:"
    echo -e "i:data input path\no:deeptools plot output path\nt:tissue"
}
while getopts ":h:i:o:t:" OPT
do
    case $OPT in
        i) data_path=${OPTARG};;
        o) result_path=${OPTARG};;
        t) tissue=${OPTARG};;        
        h) func;;
        ?) func;;
    esac
done

antibodys=(H3K27me3 H3K9me3 H3K36me3) 
ref=mm10
for antibody in ${antibodys[@]}
do 
  bam=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
  # blacklist=/storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10-blacklist.v2.bed
  samtools merge -o ${data_path}${tissue}/${antibody}/tmp.merge.bam ${bam} -@ 16
  samtools index ${data_path}${tissue}/${antibody}/tmp.merge.bam -@ 16
  window_size=1000
  gap_size=3000
  # window_size=5000
  # gap_size=10000
  e_value=100
  # bedtools bamtobed -i ${data_path}merge_bam/LLX_148_154.bam > ${data_path}bed/LLX_148_154.bed 
  mkdir ${data_path}${tissue}/${antibody}/peaks
  sicer  -t ${data_path}${tissue}/${antibody}/tmp.merge.bam  -o ${data_path}${tissue}/${antibody}/peaks  -s ${ref} -w ${window_size} -rt 16 -f 300 -egf 0.8 -fdr 0.01 -g ${gap_size} -e ${e_value} -cpu 21
  awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR, $4}' ${data_path}${tissue}/${antibody}/peaks/tmp.merge-W${window_size}-G${gap_size}.scoreisland |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/bed/${antibody}_merge-W${window_size}-G${gap_size}-E${e_value}.bed  
  # bedtools intersect -a  ${data_path}${tissue}/${antibody}/bed/${antibody}_merge-W${window_size}-G${gap_size}.bed -b ${blacklist} -v > ${data_path}${tissue}/${antibody}/bed/${antibody}_merge_rm_blacklist-W${window_size}-G${gap_size}.bed

  bw=$(ls ${data_path}${tissue}/${antibody}/bw/*nodup.bw)
  mkdir ${data_path}${tissue}/${antibody}/matrix
  computeMatrix scale-regions -S ${bw} \
    -R ${data_path}${tissue}/${antibody}/bed/${antibody}_merge-W${window_size}-G${gap_size}-E${e_value}.bed     \
    --beforeRegionStartLength 2000 --startLabel start --endLabel end \
    --regionBodyLength 2000 \
    --afterRegionStartLength 2000 \
    --numberOfProcessors 10 \
    --skipZeros -o ${data_path}${tissue}/${antibody}/matrix/${antibody}_merge_Sicer2_W${window_size}-G${gap_size}-E${e_value}.mat.gz 
  mkdir ${result_path}${tissue}
  mkdir ${result_path}${tissue}/deeptools
  plotProfile -m  ${data_path}${tissue}/${antibody}/matrix/${antibody}_merge_Sicer2_W${window_size}-G${gap_size}-E${e_value}.mat.gz   \
                -out ${result_path}${tissue}/deeptools/${antibody}_merge_Sicer2_W${window_size}-G${gap_size}-E${e_value}.pdf --startLabel start --endLabel end \
                --numPlotsPerRow 5 --perGroup --legendLocation best \
                --plotTitle "${antibody}"

  plotProfile -m  ${data_path}${tissue}/${antibody}/matrix/${antibody}_merge_Sicer2_W${window_size}-G${gap_size}-E${e_value}.mat.gz   \
                -out ${result_path}${tissue}/deeptools/${antibody}_merge_Sicer2_W${window_size}-G${gap_size}-E${e_value}_each_sample.pdf --startLabel start --endLabel end \
                --numPlotsPerRow 5  --legendLocation best \
                --plotTitle "${antibody}"

  samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
  files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
  bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}${tissue}/${antibody}/bed/${antibody}_merge-W${window_size}-G${gap_size}-E${e_value}.bed  ${data_path}${tissue}/${antibody}/bed/${antibody}_merge-W${window_size}-G${gap_size}-E${e_value}.saf
  featureCounts -p -a ${data_path}${tissue}/${antibody}/bed/${antibody}_merge-W${window_size}-G${gap_size}-E${e_value}.saf -o ${data_path}${tissue}/${antibody}/${antibody}_merge-W${window_size}-G${gap_size}-E${e_value}.counts ${files} -F SAF -T 8 
  rm ${data_path}${tissue}/${antibody}/tmp*
done
# rm ${data_path}${tissue}/${antibody}/bed/${antibody}_merge4counts.bed
# rm ${data_path}${tissue}/${antibody}/bed/${antibody}_merge.bed 
# echo ${ref}
# for sample in $samples
# do 
#     # samtools index ${data_path}${sample}.bam
#     mkdir  ${data_path}${tissue}/${antibody}/peaks/${sample}
#     sicer  -t ${data_path}${tissue}/${antibody}/bam/${sample}.nodup.bam -o ${data_path}${tissue}/${antibody}/peaks/${sample}  -s ${ref} -w 5000 -rt 1 -f 300 -egf 0.8 -fdr 0.01 -g 10000 -e 100
#     echo ${sample} done
# done
# peaks=$(ls ${data_path}${tissue}/${antibody}/peaks/*/*.scoreisland)
# cat $peaks | cut -f 1-3 |sort -k1,1 -k2,2n -|bedtools merge -i stdin | awk -v OFS='\t' 'BEGIN{num=0}{num++; print $0, "peak"num}' - |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/npc_yulab/code/bulk_atac/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/bed/${antibody}_merge4counts.bed
# bedtools intersect -a ${data_path}${tissue}/${antibody}/bed/${antibody}_merge4counts.bed -b ${blacklist} -v > ${data_path}${tissue}/${antibody}/bed/${antibody}_merge4counts_rm_blacklist.bed

# bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}${tissue}/${antibody}/bed/${antibody}_merge4counts_rm_blacklist.bed ${data_path}${tissue}/${antibody}/bed/${antibody}_merge4counts_rm_blacklist.saf

# files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
# featureCounts -p -a ${data_path}${tissue}/${antibody}/bed/${antibody}_merge4counts_rm_blacklist.saf -o ${data_path}${tissue}/${antibody}/${antibody}_merge4counts.counts ${files} -F SAF -T 8 




# macs2 callpeak -t ${data_path}merge_bam/LLX_145_151.bam -f BAMPE -n LLX_145_151 --outdir ${data_path}peaks/LLX145_151 -g mm  --keep-dup all

# cat  ${data_path}peaks/LLX145_151/LLX_145_151_peaks.narrowPeak | cut -f 1-3 |sort -k1,1 -k2,2n -|bedtools merge -i stdin | awk -v OFS='\t' 'BEGIN{num=0}{num++; print $0, "peak"num}' - |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/npc_yulab/code/bulk_atac/keep_regular_chroms.r > ${data_path}peaks/LLX145_151/LLX_145_151_peaks.narrowPeak.bed

# computeMatrix scale-regions -S ${data_path}bigWig/LLX145*.bw ${data_path}bigWig/LLX151*.bw \
#   -R  ${data_path}peaks/LLX145_151/LLX_145_151_peaks.narrowPeak.bed \
#   --beforeRegionStartLength 2000 --startLabel start --endLabel end \
#   --regionBodyLength 2000 \
#   --afterRegionStartLength 2000 \
#   --numberOfProcessors 10 \
#   --skipZeros -o ${data_path}matrix/H3K9me3_merge_macs2.mat.gz 

# plotProfile -m ${data_path}matrix/H3K9me3_merge_macs2.mat.gz  \
#               -out ${data_path}heatmap/H3K9me3_merge_LLX145_151_macs2.pdf --startLabel start --endLabel end \
#               --numPlotsPerRow 5 --perGroup --legendLocation best \
#               --plotTitle "H3K9me3_merge_LLX145_151"