ChromHMM_path=~/software/ChromHMM/
java -Xmx16G -jar ${ChromHMM_path}ChromHMM.jar CompareModels ~/projects/Aging_CUT_Tag/result/all/ChromHMM/until_ovary/25_until_ovary/emissions_25.txt \
   ~/projects/Aging_CUT_Tag/result/all/ChromHMM/until_ovary/comparedir/ compare_to_25_state_models