data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
window_size=1000
gap_size=3000
e_value=100
ref_data=~/ref_data/
tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum ileum bonemarrow)
antibodys=(H3K27me3 H3K9me3 H3K36me3)
bin_size=10kb
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        bedtools intersect -a ${ref_data}mm10_${bin_size}_bins.bed \
            -b ${data_path}${tissue}/${antibody}/bed/${antibody}_merge-W${window_size}-G${gap_size}-E${e_value}.bed -wa > ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_in-W${window_size}-G${gap_size}-E${e_value}.bed
    done
done


antibodys=(H3K27ac H3K4me1 H3K4me3)
bin_size=1kb
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        bedtools intersect -a ${ref_data}mm10_${bin_size}_bins.bed \
            -b ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_narrowpeak.bed -wa > ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_in_macs_narrowpeak.bed
    done
done

bin_size=10kb
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        bedtools intersect -a ${ref_data}mm10_${bin_size}_bins.bed \
            -b ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_narrowpeak.bed -wa > ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_in_macs_narrowpeak.bed
    done
done