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
tissue_label_change() {
    local tissue="$1"
    local tissue_label

    case "$tissue" in
        "brain")
            tissue_label="Cortex"
            ;;
        "Hip")
            tissue_label="Hippocampus"
            ;;
        "CB")
            tissue_label="Cerebellum"
            ;;
        *)
            tissue_label=$(echo "$tissue" | awk '{print toupper(substr($0, 1, 1))tolower(substr($0, 2))}')
            case "$tissue_label" in
                "Bonemarrow")
                    tissue_label="Bone Marrow"
                    ;;
                "Bat")
                    tissue_label="BAT"
                    ;;
                "Mammarygland")
                    tissue_label="Mammary Gland"
                    ;;
                "Iwat")
                    tissue_label="iWAT"
                    ;;
            esac
            ;;
    esac
    echo "$tissue_label"
}
tissue_label=$(tissue_label_change $tissue)
mkdir ${data_path}${tissue}/TEcount_subfamily/

search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/RNA_search_table.csv
cleaned_file=$(mktemp)  
cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
samples_array=$(awk -F',' -v t="$tissue_label" 'NR > 1 && ($1 == t) {print $4}' "$cleaned_file") 
IFS=$'\n' read -r -d '' -a samples < <(echo "$samples_array" && printf '\0') 

declare -A GTF_DICT  
GTF_DICT=(  
  ["hg38"]="/storage/zhangyanxiaoLab/share/gtf/hg38.gencode.v38.annotation.gtf"  
  ["hg19"]="/storage/zhangyanxiaoLab/share/gtf/hg19.gencode.v19.annotation.gtf"  
  ["mm10"]="/storage/zhangyanxiaoLab/share/gtf/mm10.gencode.vM25.annotation.gtf"  
  ["mm10_subfamily"]="/storage/zhangyanxiaoLab/share/gtf/mm10.gencode.vM25.annotation.gtf"  
  ["mm9"]="/storage/zhangyanxiaoLab/share/gtf/mm9.gencode.vM1.annotation.gtf"  
  ["panTro6"]="/storage/zhangyanxiaoLab/share/gtf/panTro6.gtf"  
  ["calJac3"]="/storage/zhangyanxiaoLab/share/gtf/calJac3.gtf"  
  ["panPan2"]="/storage/zhangyanxiaoLab/share/gtf/panPan2.gtf"  
  ["rheMac10"]="/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/rheMac10/Macaca_mulatta.Mmul_10.112.gtf"  
)  
declare -A TE_GTF_DICT  
TE_GTF_DICT=(  
  ["mm10"]="/storage/zhangyanxiaoLab/suzhuojie/ref_data/TE_reference/mm10_rmsk_TE.gtf"
  ["mm10_subfamily"]="/storage/zhangyanxiaoLab/suzhuojie/ref_data/TE_reference/mm10_rmsk_TE_subfamily.gtf"  
) 
gene_gtf=${GTF_DICT[$species]}
TE_gtf=${TE_GTF_DICT[$species]}
echo gene GTF is $gene_gtf
echo TE GTF is $TE_gtf
for sample in ${samples[@]}
do
    echo TEcount -b ${data_path}${tissue}/bam/${sample}*.sorted.bam --sortByPos --format BAM --mode multi --GTF $gene_gtf --TE $TE_gtf --project ${sample} --outdir ${data_path}${tissue}/TEcount_subfamily/
    TEcount -b ${data_path}${tissue}/bam/${sample}*.sorted.bam --sortByPos --format BAM --mode multi --verbose 3 \
        --GTF $gene_gtf --TE $TE_gtf --project ${sample} --outdir ${data_path}${tissue}/TEcount_subfamily/ &
done
wait

len=$(find ${data_path}${tissue}/TEcount_subfamily/ -type f -name "*.cntTable" | grep -v "combined.cntTable" | wc -l)  
counts=$(ls ${data_path}${tissue}/TEcount_subfamily/*.cntTable | grep -v "combined.cntTable")
paste ${counts} | cut -f 1,$(seq -s, 2 2 $((2*len))) > ${data_path}${tissue}/TEcount_subfamily/combined.cntTable