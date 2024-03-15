data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissues=(brain liver testis colon kidney lung spleen muscle pancreas)
# tissue=liver
window_size=1000
gap_size=3000
# window_size=5000
# gap_size=10000
e_value=100
for tissue in ${tissues[@]}
do 
    mkdir ${data_path}${tissue}/H3K27me3_H3K9me3_intersect
    bedtools intersect -a ${data_path}${tissue}/H3K27me3/bed/H3K27me3_merge-W${window_size}-G${gap_size}-E${e_value}.bed \
        -b ${data_path}${tissue}/H3K9me3/bed/H3K9me3_merge-W${window_size}-G${gap_size}-E${e_value}.bed -wa > ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K27me3-W${window_size}-G${gap_size}-E${e_value}.bed

    python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/H3K27me3_H3K9me3_relationship/remove_duplicate.py ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K27me3-W${window_size}-G${gap_size}-E${e_value}.bed

    sort -k 1,1 -k2,2n ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K27me3-W${window_size}-G${gap_size}-E${e_value}_unique.bed >  ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K27me3-W${window_size}-G${gap_size}-E${e_value}_unique_sort.bed

    bedtools intersect -a ${data_path}${tissue}/H3K9me3/bed/H3K9me3_merge-W${window_size}-G${gap_size}-E${e_value}.bed \
        -b ${data_path}${tissue}/H3K27me3/bed/H3K27me3_merge-W${window_size}-G${gap_size}-E${e_value}.bed -wa > ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K9me3-W${window_size}-G${gap_size}-E${e_value}.bed

    python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/H3K27me3_H3K9me3_relationship/remove_duplicate.py ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K9me3-W${window_size}-G${gap_size}-E${e_value}.bed

    sort -k 1,1 -k2,2n ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K9me3-W${window_size}-G${gap_size}-E${e_value}_unique.bed >  ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K9me3-W${window_size}-G${gap_size}-E${e_value}_unique_sort.bed

    /usr/local/bin/Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/H3K27me3_H3K9me3_relationship/H3K27me_H3K9me3_relationship_preprocess.R tissue="${tissue}"

    rm ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K27me3-W${window_size}-G${gap_size}-E${e_value}.bed \
        ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K27me3-W${window_size}-G${gap_size}-E${e_value}_unique.bed \
        ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K9me3-W${window_size}-G${gap_size}-E${e_value}.bed \
        ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K9me3-W${window_size}-G${gap_size}-E${e_value}_unique.bed

    ##### H3K27me3 or H3K9me3 increase and decrease peak intersecrt in another histone modification
    ##### Use the H3K27me3 peak as a template
    bedtools intersect -a ${data_path}${tissue}/H3K27me3/bed/H3K27me3_merge-W${window_size}-G${gap_size}-E${e_value}.bed \
        -b ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K9me3-W${window_size}-G${gap_size}-E${e_value}_up.bed -wa > ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3-W${window_size}-G${gap_size}-E${e_value}_up_in_H3K27me3.bed

    bedtools intersect -a ${data_path}${tissue}/H3K27me3/bed/H3K27me3_merge-W${window_size}-G${gap_size}-E${e_value}.bed \
        -b ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K9me3-W${window_size}-G${gap_size}-E${e_value}_down.bed -wa > ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3-W${window_size}-G${gap_size}-E${e_value}_down_in_H3K27me3.bed

    python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/H3K27me3_H3K9me3_relationship/remove_duplicate.py ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3-W${window_size}-G${gap_size}-E${e_value}_up_in_H3K27me3.bed

    python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/H3K27me3_H3K9me3_relationship/remove_duplicate.py ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3-W${window_size}-G${gap_size}-E${e_value}_down_in_H3K27me3.bed

    sort -k 1,1 -k2,2n ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3-W${window_size}-G${gap_size}-E${e_value}_up_in_H3K27me3_unique.bed >  ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3-W${window_size}-G${gap_size}-E${e_value}_up_in_H3K27me3_unique_sort.bed
    sort -k 1,1 -k2,2n ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3-W${window_size}-G${gap_size}-E${e_value}_down_in_H3K27me3_unique.bed >  ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3-W${window_size}-G${gap_size}-E${e_value}_down_in_H3K27me3_unique_sort.bed

    rm ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3-W${window_size}-G${gap_size}-E${e_value}_up_in_H3K27me3.bed \
        ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3-W${window_size}-G${gap_size}-E${e_value}_down_in_H3K27me3.bed \
        ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3-W${window_size}-G${gap_size}-E${e_value}_up_in_H3K27me3_unique.bed \
        ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K9me3-W${window_size}-G${gap_size}-E${e_value}_down_in_H3K27me3_unique.bed

    ##### Use the H3K9me3 peak as a template
    bedtools intersect -a ${data_path}${tissue}/H3K9me3/bed/H3K9me3_merge-W${window_size}-G${gap_size}-E${e_value}.bed \
        -b ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K27me3-W${window_size}-G${gap_size}-E${e_value}_up.bed -wa > ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K27me3-W${window_size}-G${gap_size}-E${e_value}_up_in_H3K9me3.bed

    bedtools intersect -a ${data_path}${tissue}/H3K9me3/bed/H3K9me3_merge-W${window_size}-G${gap_size}-E${e_value}.bed \
        -b ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/intersect_H3K27me3-W${window_size}-G${gap_size}-E${e_value}_down.bed -wa > ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K27me3-W${window_size}-G${gap_size}-E${e_value}_down_in_H3K9me3.bed

    python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/H3K27me3_H3K9me3_relationship/remove_duplicate.py ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K27me3-W${window_size}-G${gap_size}-E${e_value}_up_in_H3K9me3.bed

    python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/H3K27me3_H3K9me3_relationship/remove_duplicate.py ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K27me3-W${window_size}-G${gap_size}-E${e_value}_down_in_H3K9me3.bed

    sort -k 1,1 -k2,2n ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K27me3-W${window_size}-G${gap_size}-E${e_value}_up_in_H3K9me3_unique.bed >  ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K27me3-W${window_size}-G${gap_size}-E${e_value}_up_in_H3K9me3_unique_sort.bed
    sort -k 1,1 -k2,2n ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K27me3-W${window_size}-G${gap_size}-E${e_value}_down_in_H3K9me3_unique.bed >  ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K27me3-W${window_size}-G${gap_size}-E${e_value}_down_in_H3K9me3_unique_sort.bed

    rm ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K27me3-W${window_size}-G${gap_size}-E${e_value}_up_in_H3K9me3.bed \
        ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K27me3-W${window_size}-G${gap_size}-E${e_value}_down_in_H3K9me3.bed \
        ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K27me3-W${window_size}-G${gap_size}-E${e_value}_up_in_H3K9me3_unique.bed \
        ${data_path}${tissue}/H3K27me3_H3K9me3_intersect/H3K27me3-W${window_size}-G${gap_size}-E${e_value}_down_in_H3K9me3_unique.bed

    /usr/local/bin/Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/H3K27me3_H3K9me3_relationship/H3K27me_H3K9me3_relationship.R tissue="${tissue}"

done