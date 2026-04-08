tissue=$1
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/

files=$(ls ${data_path}${tissue}/bam/*.sorted.bam)
featureCounts -a /storage/zhangyanxiaoLab/share/gtf/mm10.gencode.vM25.annotation.gtf -p -o ${data_path}${tissue}/combined-chrM.counts ${files} -F GTF -T 10 -t exon -g gene_name 
