data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/Hippocampus_aging/
# samples=$(find ${data_path}bigWig -type f -name "*.bw" -exec basename {} \; | sed 's/.nodup.bw$//'  | sort) 
samples=(JC_R1_H3K9me3 JC_R2_H3K9me3 VC_R1_H3K9me3 VC_R2_H3K9me3)
for sample in ${samples[@]}
do 
    sf=$(grep $sample ${data_path}all_sample.qc.txt |cut -f 8)
    bamCoverage --scaleFactor $sf -b ${data_path}bam/${sample}.nodup.bam -o ${data_path}bigWig/${sample}.scaled.filt.srt.bw --outFileFormat bigwig -bs 50 --numberOfProcessors 6 --normalizeUsing RPKM &
done

for sample in ${samples[@]}
do 
    bamCoverage -b ${data_path}bam/${sample}.nodup.bam -o ${data_path}bigWig/${sample}.nodup.bw --outFileFormat bigwig -bs 50 --numberOfProcessors 6 --normalizeUsing RPKM &
done

samples=(DYQ123 DYQ124 DYQ125 DYQ126)
for sample in ${samples[@]}
do 
    sf=$(grep $sample ${data_path}all_sample.qc.txt |cut -f 8)
    bamCoverage --scaleFactor $sf -b ${data_path}bam/${sample}*.nodup.bam -o ${data_path}bigWig/${sample}.scaled.filt.srt.bw --outFileFormat bigwig --binSize 1000 --smoothLength 3000 --numberOfProcessors 6 --normalizeUsing RPKM &
done

for sample in ${samples[@]}
do 
    bamCoverage -b ${data_path}bam/${sample}*.nodup.bam -o ${data_path}bigWig/${sample}_bs1000.nodup.bw --outFileFormat bigwig --binSize 1000 --smoothLength 3000 --numberOfProcessors 6 --normalizeUsing RPKM &
done