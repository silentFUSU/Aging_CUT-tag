data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
ChromHMM_path=~/software/ChromHMM/
ref=mm10
max_jobs=2
current_jobs() {  
    jobs -rp | wc -l  
}  
mkdir -p ${result_path}all/ChromHMM/all_tissues_normal_chr/

states=(2 3 4 5 6 7 8 9 10 11 12 13 14 16 17 18 19 20 21 22 23 24 25)
for state in ${states[@]}
do
    mkdir -p ${result_path}all/ChromHMM/all_tissues_normal_chr/${state}_all_tissues
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    java -Xmx64G -jar ${ChromHMM_path}ChromHMM.jar LearnModel -b 1000 -p 16 -holdcolumnorder ${data_path}all/combined_analysis_enhancer/binarizedData_normal_chr ${result_path}all/ChromHMM/all_tissues_normal_chr/${state}_all_tissues ${state} ${ref} &
done
wait
echo all CUTTAG done