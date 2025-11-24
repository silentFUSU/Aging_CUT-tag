data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
ChromHMM_path=~/software/ChromHMM/
mkdir -p ${result_path}all/ChromHMM/all_tissues_normal_chr/comparedir/
ln -s  ${result_path}all/ChromHMM/all_tissues_normal_chr/*_all_tissues/emissions_*.txt ${result_path}all/ChromHMM/all_tissues_normal_chr/comparedir/

java -Xmx16G -jar ${ChromHMM_path}ChromHMM.jar CompareModels ${result_path}all/ChromHMM/all_tissues_normal_chr/25_all_tissues/emissions_25.txt \
   ${result_path}all/ChromHMM/all_tissues_normal_chr/comparedir/ ${result_path}all/ChromHMM/all_tissues_normal_chr/comparedir/compare_to_25_state_models
