data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
window_size=1000
gap_size=3000
e_value=100
ref_data=~/ref_data/
tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum ileum bonemarrow heart thymus stomach skin bladder tongue aorta)
bin_size=10kb
antibodys=(H3K27me3 H3K9me3 H3K36me3)
for tissue in ${tissues[@]}
do
    bash ~/projects/Aging_CUT_Tag/code/call_peak/sicer2_youngbam_callpeaks.sh ${tissue}
    for antibody in ${antibodys[@]}
    do
        bedtools intersect -a ${ref_data}mm10_${bin_size}_bins.bed \
            -b ${data_path}${tissue}/${antibody}/bed/${antibody}_young_merge-W${window_size}-G${gap_size}-E${e_value}.bed -wa > ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_in_young-W${window_size}-G${gap_size}-E${e_value}.bed
        bedtools intersect -a ${ref_data}mm10_${bin_size}_bins.bed \
            -b ${data_path}${tissue}/${antibody}/bed/${antibody}_old_only_merge-W${window_size}-G${gap_size}-E${e_value}.bed -wa > ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_in_old_only-W${window_size}-G${gap_size}-E${e_value}.bed
    
    done
done

antibodys=(H3K27ac H3K4me1 H3K4me3)
bin_size=1kb
for tissue in ${tissues[@]}
do
    bash ~/projects/Aging_CUT_Tag/code/call_peak/macs2_youngbam_callpeaks.sh ${tissue}
    for antibody in ${antibodys[@]}
    do
        bedtools intersect -a ${ref_data}mm10_${bin_size}_bins.bed \
            -b ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_young_narrowpeak.bed -wa > ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_in_young_macs_narrowpeak.bed
                bedtools intersect -a ${ref_data}mm10_${bin_size}_bins.bed \
            -b ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_old_only_narrowpeak.bed -wa > ${data_path}${tissue}/${antibody}/bed/${antibody}_${bin_size}_in_old_only_macs_narrowpeak.bed
    done
done