tissue=$1
antibodys=(H3K27me3)
for antibody in ${antibodys[@]}
do
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "domain"NR}' ${data_path}${tissue}/${antibody}/peaks/edd/edd_peaks_fdr05.bed > ${data_path}${tissue}/${antibody}/peaks/edd/age_domain.bed
    files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
    bed=${data_path}${tissue}/${antibody}/peaks/edd/age_domain.bed
    saf=${data_path}${tissue}/${antibody}/peaks/edd/age_domain.saf
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${bed} ${saf}
    featureCounts -p -C -a ${saf} -o ${data_path}${tissue}/${antibody}/${antibody}_edd_peaks_fdr05.counts ${files} -F SAF -T 8 
done