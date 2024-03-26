data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
ref=mm10  
tissue=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum bonemarrow)
antibodys=(H3K27me3 H3K9me3 H3K36me3)
# antibodys=(H3K4me1)
for antibody in ${antibodys[@]}
do 
    bam=()  
    for t in ${tissue[@]}
    do  
        bam+=($(ls ${data_path}${t}/${antibody}/bam/*.bam))  
    done  
    featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_10kb_bins.saf -o ${data_path}all/${antibody}/merge-10kb_bins.counts ${bam[@]} -F SAF -T 8 
done

antibodys=(H3K4me1 H3K4me3 H3K27ac)
for antibody in ${antibodys[@]}
do 
    bam=()  
    for t in ${tissue[@]}
    do  
        bam+=($(ls ${data_path}${t}/${antibody}/bam/*.bam))  
    done  
    featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_1kb_bins.saf -o ${data_path}all/${antibody}/merge-1kb_bins.counts ${bam[@]} -F SAF -T 8 
done