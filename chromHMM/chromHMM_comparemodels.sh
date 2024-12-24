data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
ChromHMM_path=~/software/ChromHMM/
mkdir -p ${result_path}all/ChromHMM/all_tissues/comparedir/
ln -s  ${result_path}all/ChromHMM/all_tissues/*_all_tissues/emissions_*.txt ${result_path}all/ChromHMM/all_tissues/comparedir/

java -Xmx16G -jar ${ChromHMM_path}ChromHMM.jar CompareModels ${result_path}all/ChromHMM/all_tissues/25_all_tissues/emissions_25.txt \
   ${result_path}all/ChromHMM/all_tissues/comparedir/ ${result_path}all/ChromHMM/all_tissues/comparedir/compare_to_25_state_models
