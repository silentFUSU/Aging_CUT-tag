data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
antibodys=(H3K27me3)
ref=mm10
tissues=$1
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
        featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS_10kb.saf -o ${data_path}${tissue}/${antibody}/${antibody}_gene_TSS_10kb.counts ${files} -F SAF -T 8 
    done
done