data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
ChromHMM_path=~/software/ChromHMM/
ref=mm10
mkdir -p ${result_path}all/ChromHMM/all_tissues

states=(2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25)
for state in ${states[@]}
do
    mkdir -p ${result_path}all/ChromHMM/all_tissues/${state}_all_tissues
    java -Xmx64G -jar ${ChromHMM_path}ChromHMM.jar LearnModel -b 1000 -p 16 -holdcolumnorder ${data_path}all/combined_analysis_enhancer/binarizedData ${result_path}all/ChromHMM/all_tissues/${state}_all_tissues ${state} ${ref}
done