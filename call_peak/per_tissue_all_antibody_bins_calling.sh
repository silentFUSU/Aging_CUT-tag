tissues=(brain liver testis colon kidney lung spleen muscle Hip cecum bonemarrow heart thymus stomach skin aorta tongue bladder CB jejunum uterus ovary ileum pancreas BAT mammarygland)
CUTTAG_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
ATAC_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
for tissue in ${tissues[@]}
do
    CUTTAG_files=$(ls ${CUTTAG_path}${tissue}/{H3K27me3,H3K9me3,H3K36me3,H3K27ac,H3K4me3,H3K4me1}/bam/*.bam)
    ATAC_files=$(ls ${ATAC_path}${tissue}/ATAC/bam/*.bam)
    CUTTAG_files_array=($CUTTAG_files)  
    ATAC_files_array=($ATAC_files)  
    files=("${CUTTAG_files_array[@]}" "${ATAC_files_array[@]}")  
    echo "featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_10kb_bins.saf -o ${CUTTAG_path}${tissue}/all_antibodys_10kb_bins.counts ${files[@]} -F SAF -T 8 "
    featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_10kb_bins.saf -o ${CUTTAG_path}${tissue}/all_antibodys_10kb_bins.counts ${files[@]} -F SAF -T 8 
done

