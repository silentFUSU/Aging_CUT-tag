deep_data_path=/storage/zhangyanxiaoLab/zhangyanxiao/software/obsutil_linux_amd64_5.2.10/tmp/YuLab/20241016/
shallow_data_path=/storage/zhangyanxiaoLab/fastq/2024/2024-09-25-Jiangbei-YuLab/
target_data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20241016_WJH_HiC/
samples=(WJH-103-Lung WJH-109-Lung)
for sample in ${samples[@]}
do
    mkdir ${target_data_path}fastq/${sample}
    f1=${deep_data_path}*${sample}*/*${sample}*_1.fq.gz
    f2=${shallow_data_path}*${sample}*/${sample}*_R1_001.fastq.gz
    f3=${target_data_path}fastq/${sample}/${sample}_R1.fastq.gz
    echo 'zcat' ${f1} ${f2} '| gzip -> ' ${f3}
    zcat ${f1} ${f2} | gzip -> ${f3} &
    f1=${deep_data_path}*${sample}*/*${sample}*_2.fq.gz
    f2=${shallow_data_path}*${sample}*/${sample}*_R2_001.fastq.gz
    f3=${target_data_path}fastq/${sample}/${sample}_R2.fastq.gz
    echo 'zcat' ${f1} ${f2} '| gzip -> ' ${f3} 
    zcat ${f1} ${f2} | gzip -> ${f3} &
done
wait
echo merge done