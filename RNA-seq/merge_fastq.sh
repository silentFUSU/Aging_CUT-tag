deep_path=/storage/zhangyanxiaoLab/fastq/2024/2024-05-20-Jiangbei-SZJ/
shallow_path=/storage/zhangyanxiaoLab/fastq/2024/2024-05-10-Jiangbei-SZJ/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240520_SZJ/

mkdir ${result_path}
mkdir ${result_path}rawdata/
samples=$(find ${deep_path}  -type d -name "SZJ*" -exec basename {} \;) 
# samples=LLX494
for sample in $samples
do
    f1=${deep_path}${sample}/${sample}_*_R1*
    f2=${shallow_path}${sample}/${sample}_*_R1*
    f3=${result_path}rawdata/${sample}_R1.fastq.gz
    zcat ${f1} ${f2} | gzip -> ${f3} &
    echo ${f1} ${f2} ${f3} combine
    f1=${deep_path}${sample}/${sample}_*_R2*
    f2=${shallow_path}${sample}/${sample}_*_R2*
    f3=${result_path}rawdata/${sample}_R2.fastq.gz
    zcat ${f1} ${f2} | gzip -> ${f3} &
    echo ${f1} ${f2} ${f3} combine
done
wait
echo all done