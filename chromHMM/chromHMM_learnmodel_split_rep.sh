data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
ChromHMM_path=~/software/ChromHMM/
ref=mm10
reps=(rep1 rep2)
state=15
for rep in ${reps[@]}
do
    mkdir -p ${result_path}all/ChromHMM/all_tissues_normal_chr_split_rep/${state}_all_tissues_${rep}/
    java -Xmx64G -jar ${ChromHMM_path}ChromHMM.jar LearnModel -b 1000 -p 16 -holdcolumnorder -init load -m /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/ChromHMM/all_tissues_normal_chr/15_all_tissues/model_15.txt -r 50 ${data_path}all/combined_analysis_enhancer/binarizedData_normal_chr_split_rep/${rep} ${result_path}all/ChromHMM/all_tissues_normal_chr_split_rep/${state}_all_tissues_${rep}/ ${state} ${ref} &
done