data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
ChromHMM_path=~/software/ChromHMM/
ref=mm10
java -Xmx64G -jar ${ChromHMM_path}ChromHMM.jar BinarizeBed -b 1000 ${ChromHMM_path}CHROMSIZES/${ref}.txt \
    ${data_path}all/chromHMM_with_ATAC/chromHMMbed/ \
    ${data_path}all/chromHMM_with_ATAC/cellmarkfiletable.txt ${data_path}all/chromHMM_with_ATAC/binarizedData/ &
wait
echo all done