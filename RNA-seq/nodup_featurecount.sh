tissue=$1
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
samples=$(find ${data_path}${tissue}/bam/ -type f -name "*.nodup.bam" -exec basename {} \; | sed 's/\.nodup.bam//' | sort)  
for sample in ${samples[@]}
do
    featureCounts -a /storage/zhangyanxiaoLab/share/gtf/mm10.gencode.vM25.annotation.gtf -p -o ${data_path}${tissue}/counts/${sample}.nodup.counts  ${data_path}${tissue}/bam/${sample}*.nodup.bam -F GTF -T 10 -t exon -g gene_name &
done


len=$(find ${data_path}${tissue}/counts/ -type f -name "*.nodup.counts" | wc -l)  
counts=$(ls ${data_path}${tissue}/counts/*.nodup.counts)
paste ${counts} | cut -f 1-6,$(seq -s, 7 7 $((7*len))) | grep -v 'chrM' > ${data_path}${tissue}/combined-chrM.nodup.counts