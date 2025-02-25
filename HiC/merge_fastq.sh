deep_data_path=/storage/zhangyanxiaoLab/zhangyanxiao/software/obsutil_linux_amd64_5.2.10/tmp/LLX/20250220-1/
shallow_data_path=/storage/zhangyanxiaoLab/fastq/2025/2025-01-17-Jiangbei-LLX/
target_data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data_transposon2/20250220_LLX_CUTTag/
# samples=(WJH_stomach_96 WJH_stomach_99 WJH_stomach_108 WJH_stomach_112 WJH_thymus_230 WJH_thymus_240 WJH_thymus_241 WJH_thymus_245)
samples=(WJH_stomach_99)
mkdir -p ${target_data_path}
mkdir -p ${target_data_path}fastq/

# for sample in ${samples[@]}
# do
#     mkdir -p ${target_data_path}fastq/
#     mkdir -p ${target_data_path}fastq/${sample}
#     f1=${deep_data_path}*${sample}*/*${sample}*_1.fq.gz
#     f2=${shallow_data_path}*${sample}*/*${sample}*_1.fq.gz
#     f3=${target_data_path}fastq/${sample}/${sample}_R1.fastq.gz
#     echo 'cat' ${f1} ${f2} ' > ' ${f3}
#     cat ${f1} ${f2} > ${f3} &
#     f1=${deep_data_path}*${sample}*/*${sample}*_2.fq.gz
#     f2=${shallow_data_path}*${sample}*/*${sample}*_2.fq.gz
#     f3=${target_data_path}fastq/${sample}/${sample}_R2.fastq.gz
#     echo 'cat' ${f1} ${f2} ' > ' ${f3} 
#     cat ${f1} ${f2} > ${f3} &
# done

for sample in ${samples[@]}
do
    mkdir -p ${target_data_path}fastq/
    mkdir -p ${target_data_path}fastq/${sample}
    f1=${deep_data_path}*${sample}*R1.raw.fastq.gz
    f2=${shallow_data_path}*${sample}*raw.1.fastq.gz
    f3=${target_data_path}fastq/${sample}/${sample}_R1.fastq.gz
    echo 'cat' ${f1} ${f2} ' > ' ${f3}
    cat ${f1} ${f2} > ${f3} &
    f1=${deep_data_path}*${sample}*R2.raw.fastq.gz
    f2=${shallow_data_path}*${sample}*raw.2.fastq.gz
    f3=${target_data_path}fastq/${sample}/${sample}_R2.fastq.gz
    echo 'cat' ${f1} ${f2} ' > ' ${f3} 
    cat ${f1} ${f2} > ${f3} &
done
wait
echo merge done