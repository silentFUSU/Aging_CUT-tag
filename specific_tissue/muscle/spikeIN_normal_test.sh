ln -s /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/muscle/H3K36me3/bam/* /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/muscle_cut_tag_test/H3K36me3/bam/
ln -s /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/muscle/H3K4me3/bam/* /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/muscle_cut_tag_test/H3K4me3/bam/
ln -s /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240813_CKJ_CUTTAG/bam/CKJ{060,061,062,063}* /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/muscle_cut_tag_test/H3K36me3/bam/
ln -s /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240813_CKJ_CUTTAG/bam/CKJ{064,065,066,067}* /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/muscle_cut_tag_test/H3K4me3/bam/
files=$(ls /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/muscle_cut_tag_test/H3K36me3/bam/*.bam)
featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_10kb_bins.saf -o /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/muscle_cut_tag_test/H3K36me3/10kb_bins.counts ${files} -F SAF -T 8 
files=$(ls /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/muscle_cut_tag_test/H3K4me3/bam/*.bam)
featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_1kb_bins.saf -o /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/muscle_cut_tag_test/H3K4me3/1kb_bins.counts ${files} -F SAF -T 8 

data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/muscle_cut_tag_test/
antibody=H3K36me3
mkdir ${data_path}${antibody}/peaks ${data_path}${antibody}/featurecount/ ${data_path}${antibody}/bed/
window_size=1000
gap_size=3000
e_value=100
ref=mm10
blacklist=~/ref_data/mm10-blacklist.v2.bed
samples=(CKJ060 CKJ061 CKJ062 CKJ063 LLX280 LLX286 LLX292 LLX298)
for sample in ${samples[@]}
do
    sicer  -t ${data_path}${antibody}/bam/${sample}*.bam -o ${data_path}${antibody}/peaks  -s ${ref} -w ${window_size} -rt 16 -f 300 -egf 0.8 -fdr 0.01 -g ${gap_size} -e ${e_value} -cpu 21
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR, $4}' ${data_path}${antibody}/peaks/${sample}.nodup-W${window_size}-G${gap_size}.scoreisland |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${antibody}/bed/${sample}.nodup-W${window_size}-G${gap_size}.bed 

    bedtools intersect -a ${data_path}${antibody}/bed/${sample}.nodup-W${window_size}-G${gap_size}.bed \
        -b ${blacklist} -v > ${data_path}${antibody}/bed/${sample}_remove_blacklist-W${window_size}-G${gap_size}.bed
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}${antibody}/bed/${sample}_remove_blacklist-W${window_size}-G${gap_size}.bed  ${data_path}${antibody}/bed/${sample}_remove_blacklist-W${window_size}-G${gap_size}.saf 
    featureCounts -p -a ${data_path}${antibody}/bed/${sample}_remove_blacklist-W${window_size}-G${gap_size}.saf  -o ${data_path}${antibody}/featurecount/${sample}_remove_blacklist.counts ${data_path}${antibody}/bam/${sample}*.bam  -F SAF -T 8 &
done

data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/muscle_cut_tag_test/
antibody=H3K4me3
mkdir ${data_path}${antibody}/peaks ${data_path}${antibody}/featurecount/ ${data_path}${antibody}/bed/ ${data_path}/${antibody}/peaks/macs_narrowpeak
ref=mm
blacklist=~/ref_data/mm10-blacklist.v2.bed
samples=(CKJ065 CKJ066 CKJ067 LLX284 LLX290 LLX296 LLX302)
for sample in ${samples[@]}
do
    macs2 callpeak -t ${data_path}${antibody}/bam/${sample}*.bam  -f BAMPE -n ${sample} --outdir ${data_path}/${antibody}/peaks/macs_narrowpeak -g ${ref} --nomodel -q 0.0001  --keep-dup all
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' ${data_path}/${antibody}/peaks/macs_narrowpeak/${sample}_peaks.narrowPeak |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}/${antibody}/bed/${sample}_peaks.bed 
    
    bedtools intersect -a ${data_path}/${antibody}/bed/${sample}_peaks.bed \
        -b ${blacklist} -v > ${data_path}/${antibody}/bed/${sample}_remove_blacklist_peaks.bed
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}/${antibody}/bed/${sample}_remove_blacklist_peaks.bed  ${data_path}/${antibody}/bed/${sample}_remove_blacklist_peaks.saf 
    featureCounts -p -a  ${data_path}/${antibody}/bed/${sample}_remove_blacklist_peaks.saf -o ${data_path}${antibody}/featurecount/${sample}_remove_blacklist.counts ${data_path}${antibody}/bam/${sample}*.bam  -F SAF -T 8 &
done