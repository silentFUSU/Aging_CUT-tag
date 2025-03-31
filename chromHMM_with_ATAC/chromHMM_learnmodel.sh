data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
ChromHMM_path=~/software/ChromHMM/
ref=mm10
max_jobs=5
current_jobs() {  
    jobs -rp | wc -l  
}  
mkdir -p ${result_path}all/ChromHMM_with_ATAC/
mkdir -p ${result_path}all/ChromHMM_with_ATAC/all_tissues/

states=(10 11 12 13 14 15 16 17 18 19 20)
for state in ${states[@]}
do
    mkdir -p ${result_path}all/ChromHMM_with_ATAC/all_tissues/${state}_all_tissues
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    java -Xmx64G -jar ${ChromHMM_path}ChromHMM.jar LearnModel -b 1000 -p 16 -holdcolumnorder ${data_path}all/chromHMM_with_ATAC/binarizedData ${result_path}all/ChromHMM_with_ATAC/all_tissues/${state}_all_tissues ${state} ${ref} &
done
wait
echo all CUTTAG done