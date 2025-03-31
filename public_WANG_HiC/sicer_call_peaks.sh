data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/cellular_aging_GSE133292/Chip_seq/
samples=(SRR13274560 SRR13274561 SRR13274566 SRR13274567)
ref=hg38
antibody=H3K9me3
bam=($(ls ${data_path}bam/{SRR13274560,SRR13274561,SRR13274566,SRR13274567}.nodup.bam))
declare -a young
declare -a old
young=()
old=()
young+=("${bam[0]}")
young+=("${bam[1]}")
old+=("${bam[2]}")
old+=("${bam[3]}")
echo ${young[@]}
echo ${old[@]}
echo "samtools merge -o ${data_path}tmp.young.merge.bam ${young[@]} -@ 16"
echo "samtools merge -o ${data_path}tmp.old.merge.bam ${old[@]} -@ 16"
samtools merge -f -o ${data_path}tmp.young.merge.bam ${young[@]} -@ 16 &
samtools merge -f -o ${data_path}tmp.old.merge.bam ${old[@]} -@ 16 &
wait
samtools index ${data_path}tmp.young.merge.bam -@ 16 &
samtools index ${data_path}tmp.old.merge.bam -@ 16 &
wait
window_size=5000
gap_size=10000
e_value=100
data_path}peaks/
sicer  -t ${data_path}tmp.young.merge.bam  -o ${data_path}peaks  -s ${ref} -w ${window_size} -rt 16 -f 300 -egf 0.8 -fdr 0.01 -g ${gap_size} -e ${e_value} -cpu 21 &
sicer  -t ${data_path}tmp.old.merge.bam  -o ${data_path}peaks  -s ${ref} -w ${window_size} -rt 16 -f 300 -egf 0.8 -fdr 0.01 -g ${gap_size} -e ${e_value} -cpu 21 &
wait
awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR, $4}' ${data_path}peaks/tmp.young.merge-W${window_size}-G${gap_size}.scoreisland |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}bed/${antibody}_young_merge-W${window_size}-G${gap_size}-E${e_value}.bed  
awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR, $4}' ${data_path}peaks/tmp.old.merge-W${window_size}-G${gap_size}.scoreisland |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}bed/${antibody}_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed  
rm ${data_path}tmp*
bedtools subtract -a ${data_path}bed/${antibody}_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed -b ${data_path}bed/${antibody}_young_merge-W${window_size}-G${gap_size}-E${e_value}.bed > ${data_path}bed/${antibody}_old_only_merge-W${window_size}-G${gap_size}-E${e_value}.bed
bedtools intersect -a ${data_path}bed/${antibody}_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed -b ${data_path}bed/${antibody}_young_merge-W${window_size}-G${gap_size}-E${e_value}.bed >  ${data_path}bed/${antibody}_young_old_intersect-W${window_size}-G${gap_size}-E${e_value}.bed
cat  ${data_path}bed/${antibody}_young_merge-W${window_size}-G${gap_size}-E${e_value}.bed ${data_path}bed/${antibody}_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge  > ${data_path}bed/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed

bin_size=10kb
bedtools intersect -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/hg38/hg38_10kb_bins.bed \
            -b ${data_path}bed/${antibody}_young_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed -wa > ${data_path}bed/${antibody}_${bin_size}_in_young_old_merge-W${window_size}-G${gap_size}-E${e_value}.bed