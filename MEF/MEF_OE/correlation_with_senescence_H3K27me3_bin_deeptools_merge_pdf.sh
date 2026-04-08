conditions=(MEF_Bmi1 MEF_Cbx2 MEF_Cbx7)
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/MEF_OE/
quadrants=(first second third fourth)
for condition in ${conditions[@]}
do
    files=()
    for quadrant in ${quadrants[@]}
    do
        files+=("${result_path}plot_bin/H3K27me3_${condition}_${quadrant}.pdf")
        files+=("${result_path}plot_bin/H2AK119ub1_${condition}_${quadrant}.pdf")
        files+=("${result_path}plot_bin/senescence_H3K27me3_in_${condition}_${quadrant}.pdf")
    done
    /storage/zhangyanxiaoLab/suzhuojie/software/pdfjam-3.11/bin/pdfjam ${files[@]} --nup 3x4 --landscape --outfile ${result_path}/plot_bin/${condition}.pdf
done
