# bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/RNA-seq/TEcount.sh -i /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/ -g mm10 -t pancreas
func() {
    echo "Usage:"
    echo -e "i:data input path\ng:species\nt:tissue"
}
while getopts ":h:i:g:t:" OPT
do
    case $OPT in
        i) data_path=${OPTARG};;
        g) species=${OPTARG};;
        t) tissue=${OPTARG};;        
        h) func;;
        ?) func;;
    esac
done
# data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
# tissue=pancreas
# species=mm10
mkdir ${data_path}${tissue}/TEcount/
samples=($(find ${data_path}${tissue}/bam -type f -name "*.sorted.bam" | grep -E '/((LLX|CKJ|SRR)[0-9]+)' | awk -F'/' '{split($NF,a,"[_|.]"); print a[1]}' ))
declare -A GTF_DICT  
GTF_DICT=(  
  ["hg38"]="/storage/zhangyanxiaoLab/share/gtf/hg38.gencode.v38.annotation.gtf"  
  ["hg19"]="/storage/zhangyanxiaoLab/share/gtf/hg19.gencode.v19.annotation.gtf"  
  ["mm10"]="/storage/zhangyanxiaoLab/share/gtf/mm10.gencode.vM25.annotation.gtf"  
  ["mm9"]="/storage/zhangyanxiaoLab/share/gtf/mm9.gencode.vM1.annotation.gtf"  
  ["panTro6"]="/storage/zhangyanxiaoLab/share/gtf/panTro6.gtf"  
  ["calJac3"]="/storage/zhangyanxiaoLab/share/gtf/calJac3.gtf"  
  ["panPan2"]="/storage/zhangyanxiaoLab/share/gtf/panPan2.gtf"  
  ["rheMac10"]="/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/rheMac10/Macaca_mulatta.Mmul_10.112.gtf"  
)  
declare -A TE_GTF_DICT  
TE_GTF_DICT=(  
  ["mm10"]="/storage/zhangyanxiaoLab/suzhuojie/ref_data/TE_reference/mm10_rmsk_TE.gtf"  
) 
gene_gtf=${GTF_DICT[$species]}
TE_gtf=${TE_GTF_DICT[$species]}
echo gene GTF is $gene_gtf
echo TE GTF is $TE_gtf
for sample in ${samples[@]}
do
    echo TEcount -b ${data_path}${tissue}/bam/${sample}*.sorted.bam --sortByPos --format BAM --mode multi --GTF $gene_gtf --TE $TE_gtf --project ${sample} --outdir ${data_path}${tissue}/TEcount/
    TEcount -b ${data_path}${tissue}/bam/${sample}*.sorted.bam --sortByPos --format BAM --mode multi \
        --GTF $gene_gtf --TE $TE_gtf --project ${sample} --outdir ${data_path}${tissue}/TEcount/ &
done
wait

len=$(find ${data_path}${tissue}/TEcount/ -type f -name "*.cntTable" | grep -v "combined.cntTable" | wc -l)  
counts=$(ls ${data_path}${tissue}/TEcount/*.cntTable | grep -v "combined.cntTable")
paste ${counts} | cut -f 1,$(seq -s, 2 2 $((2*len))) > ${data_path}${tissue}/TEcount/combined.cntTable