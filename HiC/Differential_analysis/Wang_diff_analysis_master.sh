tissue=$1
code_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/HiC/
bash ${code_path}build_HiCpro_matrix.sh ${tissue} 200000
bash ${code_path}sparse2dense.sh ${tissue} 200000
python ${code_path}Differential_analysis/step1_diff_interaction_for_intrachrom.py -t ${tissue} -r 200000
python ${code_path}Differential_analysis/step2_diff_interaction_for_intrachrom.py -t ${tissue} -r 200000
/usr/local/lib64/R/bin/Rscript ${code_path}Differential_analysis/step3_diff_interaction_for_intrachrom.R -t ${tissue} -r 200000 &
bash ${code_path}with_histone_modification/interaction_change_in_histone_condition.sh $tissue &
wait
echo preprocess done