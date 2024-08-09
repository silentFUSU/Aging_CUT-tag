tissue=$1
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
samples=$(find ${data_path}${tissue}/bam/ -type f -name "*.nodup.bam" -exec basename {} \; | sed 's/\.nodup.bam//' | sort)  
for sample in ${samples[@]}
do
    bamCoverage -b ${data_path}${tissue}/bam/${sample}*.nodup.bam  -o ${data_path}${tissue}/bw/${sample}.nodup.bw --outFileFormat bigwig -bs 50 --numberOfProcessors 10 --normalizeUsing RPKM &
done   
