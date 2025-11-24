data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
antibodys=(H3K27me3 H3K9me3 H3K36me3 H3K27ac H3K4me1 H3K4me3)
# tissues=(brain liver testis colon kidney lung spleen muscle pancreas)
# tissues=(testis colon kidney lung spleen pancreas muscle Hip)
tissues=$1
ref=mm10
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
        featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_10kb_bins.saf -o ${data_path}${tissue}/${antibody}/${antibody}_10kb_bins.counts ${files} -F SAF -T 8 
    done
done

antibodys=(H3K27ac H3K4me1 H3K4me3)
# tissues=(bonemarrow)
# antibodys=(ATAC)
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
        featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_1kb_bins.saf -o ${data_path}${tissue}/${antibody}/${antibody}_1kb_bins.counts ${files} -F SAF -T 8 
    done
done

data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
antibodys=(ATAC)
# tissues=(bonemarrow)
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
        featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_10kb_bins.saf -o ${data_path}${tissue}/${antibody}/${antibody}_10kb_bins.counts ${files} -F SAF -T 8 
    done
done

# tissues=(bonemarrow)
# antibodys=(ATAC)
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
        featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_1kb_bins.saf -o ${data_path}${tissue}/${antibody}/${antibody}_1kb_bins.counts ${files} -F SAF -T 8 
    done
done