tissue=cecum
# antibodys=(H3K27me3 H3K36me3 H3K27ac H3K9me3 H3K4me3 H3K4me1 ATAC)
# antibodys=(H3K27ac H3K4me1)
# antibodys=(H3K27me3 H3K36me3 H3K27ac H3K9me3 H3K4me3 H3K4me1)
antibodys=(ATAC)
dirs=(bam bed bw)
# data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
mkdir ${data_path}${tissue}
for antibody in ${antibodys[@]}
do 
    mkdir ${data_path}${tissue}/${antibody}/
    for dir in ${dirs[@]}
    do
        mkdir ${data_path}${tissue}/${antibody}/${dir}
    done 
done

