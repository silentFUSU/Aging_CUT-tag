tissues=(lung brain liver CB)
code_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/HiC/with_histone_modification/
age=24m
for tissue in ${tissues[@]}
do
    bash ${code_path}interaction_change_in_histone_condition.sh ${tissue} ${age}
done