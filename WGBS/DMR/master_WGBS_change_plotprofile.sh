tissues=(mammarygland)
code_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/WGBS/DMR/
for tissue in ${tissues[@]}
do
    bash ${code_path}WGBS_change_in_histone_change_plotprofile.sh ${tissue}
    bash ${code_path}WGBS_change_in_DEG_TSS_plotprofile.sh ${tissue}
    bash ${code_path}WGBS_change_in_DEG_body_plotprofile.sh ${tissue}
done