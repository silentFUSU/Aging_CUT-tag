deep_path=/storage/zhangyanxiaoLab/zhangyanxiao/software/oss_download/download/240805-103234/Data/
shallow_path=/storage/zhangyanxiaoLab/zhangyanxiao/software/oss_download/download/240805-103234/Data/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240730_LLX637_638_RNA/

mkdir ${result_path}
mkdir ${result_path}rawdata/
# samples=$(find ${deep_path}  -type d -name "SZJ*" -exec basename {} \;) 
# # samples=LLX494
# for sample in $samples
# do
    f1=${deep_path}LLX637/LLX637*_R1*.gz
    f2=${shallow_path}LLX638/LLX638*_R1*.gz
    f3=${result_path}rawdata/LLX637_R1.fastq.gz
    zcat ${f1} ${f2} | gzip -> ${f3} &
    echo ${f1} ${f2} ${f3} combine
    f1=${deep_path}LLX637/LLX637*_R2*.gz
    f2=${shallow_path}LLX638/LLX638*_R2*.gz
    f3=${result_path}rawdata/LLX637_R2.fastq.gz
    zcat ${f1} ${f2} | gzip -> ${f3} &
    echo ${f1} ${f2} ${f3} combine
# done
wait
echo all done