deep_data_path=/storage/zhangyanxiaoLab/fastq/2024/2024-11-11-Jiangbei-YuLab/20241127/
shallow_data_path=/storage/zhangyanxiaoLab/fastq/2024/2024-11-05-Jiangbei-YuLab/
target_data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data_transposon2/20241220_CB_rerun_HiC/
samples=(WJH-Cerebellum-99 WJH-Cerebellum-108)
mkdir -p ${target_data_path}
mkdir -p ${target_data_path}fastq/

for sample in ${samples[@]}
do
    mkdir -p ${target_data_path}fastq/
    mkdir -p ${target_data_path}fastq/${sample}
    f1=${deep_data_path}*${sample}*/*${sample}*_1.fq.gz
    f2=${shallow_data_path}*${sample}*/*${sample}*_1.fq.gz
    f3=${target_data_path}fastq/${sample}/${sample}_R1.fastq.gz
    echo 'cat' ${f1} ${f2} ' > ' ${f3}
    cat ${f1} ${f2} > ${f3} &
    f1=${deep_data_path}*${sample}*/*${sample}*_2.fq.gz
    f2=${shallow_data_path}*${sample}*/*${sample}*_2.fq.gz
    f3=${target_data_path}fastq/${sample}/${sample}_R2.fastq.gz
    echo 'cat' ${f1} ${f2} ' > ' ${f3} 
    cat ${f1} ${f2} > ${f3} &
done
wait
echo merge done