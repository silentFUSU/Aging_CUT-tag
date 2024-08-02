data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
tissues=(tongue)
antibodys=(H3K27me3 H3K36me3 H3K27ac H3K9me3 H3K4me3 H3K4me1)
for tissue in ${tissues[@]}
do
    mkdir ~/public_html/Aging_CUT_Tag/${tissue}
    for antibody in ${antibodys[@]}
    do 
        mkdir ~/public_html/Aging_CUT_Tag/${tissue}/${antibody}/
        mkdir ~/public_html/Aging_CUT_Tag/${tissue}/${antibody}/bw
        mkdir ~/public_html/Aging_CUT_Tag/${tissue}/${antibody}/bed
    done
    for antibody in ${antibodys[@]}
    do
    ln -s ${data_path}${tissue}/${antibody}/bw/* ~/public_html/Aging_CUT_Tag/${tissue}/${antibody}/bw
    ln -s ${data_path}${tissue}/${antibody}/bed/* ~/public_html/Aging_CUT_Tag/${tissue}/${antibody}/bed
    done
done