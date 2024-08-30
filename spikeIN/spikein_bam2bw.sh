data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240827_CKJ_CUTTAG/
# samples=$(find ${data_path}bigWig -type f -name "*.bw" -exec basename {} \; | sed 's/.nodup.bw$//'  | sort) 
samples=(CKJ074 CKJ075 CKJ076 CKJ077 CKJ078 CKJ079)
for sample in ${samples[@]}
do 
    sf=$(grep $sample ${data_path}all_sample.qc.txt |cut -f 8)
    bamCoverage --scaleFactor $sf -b ${data_path}bam/${sample}.nodup.bam -o ${data_path}bigWig/${sample}.scaled.filt.srt.bw --outFileFormat bigwig -bs 50 --numberOfProcessors 6 --normalizeUsing RPKM &
done
for sample in ${samples[@]}
do 
    bamCoverage -b ${data_path}bam/${sample}.nodup.bam -o ${data_path}bigWig/${sample}.nodup.bw --outFileFormat bigwig -bs 50 --numberOfProcessors 6 --normalizeUsing RPKM &
done

samples=(CKJ068 CKJ069 CKJ070 CKJ071 CKJ072 CKJ073)
for sample in ${samples[@]}
do 
    sf=$(grep $sample ${data_path}all_sample.qc.txt |cut -f 8)
    bamCoverage --scaleFactor $sf -b ${data_path}bam/${sample}.nodup.bam -o ${data_path}bigWig/${sample}.scaled.filt.srt.bw --outFileFormat bigwig --binSize 1000 --smoothLength 3000 --numberOfProcessors 6 --normalizeUsing RPKM &
done

for sample in ${samples[@]}
do 
    bamCoverage -b ${data_path}bam/${sample}.nodup.bam -o ${data_path}bigWig/${sample}.nodup.bw --outFileFormat bigwig --binSize 1000 --smoothLength 3000 --numberOfProcessors 6 --normalizeUsing RPKM &
done