result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/RNA/TE/TE_with_histone/
antibodys=(H3K9me3 H3K27me3 H3K36me3 H3K27ac H3K4me3 H3K4me1)
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
class=ERV1_kmeans1
software_path=/storage/zhangyanxiaoLab/suzhuojie/software/pdfjam-3.11/bin/
mkdir -p ${result_path}per_tissue
for tissue in ${tissues[@]}
do
    pdfs=()
    for antibody in ${antibodys[@]}
    do
        pdfs+=(${result_path}${antibody}/plot/${tissue}_${antibody}_change_in_${class}_regions.pdf)
    done
    ${software_path}pdfjam --nup 3x2 ${pdfs[@]} -o ${result_path}per_tissue/${tissue}_change_in_${class}_regions.pdf
    convert -density 150  ${result_path}per_tissue/${tissue}_change_in_${class}_regions.pdf ${result_path}per_tissue/${tissue}_change_in_${class}_regions.png
    rm ${result_path}per_tissue/${tissue}_change_in_${class}_regions.pdf
done