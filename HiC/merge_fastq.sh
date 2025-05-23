deep_data_path=/storage/zhangyanxiaoLab/huangmin/software1/obsutil_linux_amd64_5.2.10/tmp/SZJ/20250519/
shallow_data_path=/storage/zhangyanxiaoLab/fastq/2025/2025-04-15-Jiangbei-DYQ/
target_data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data_transposon2/20250519_DYQ_WGBS/
# samples=(WJH_stomach_96 WJH_stomach_99 WJH_stomach_108 WJH_stomach_112 WJH_thymus_230 WJH_thymus_240 WJH_thymus_241 WJH_thymus_245)
samples=(DYQ171 DYQ172)
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
    f1=${deep_data_path}${sample}/*${sample}*R1*fastq.gz
    f2=${shallow_data_path}${sample}/*${sample}*R1*fastq.gz
    f3=${target_data_path}fastq/${sample}/${sample}_R1.fastq.gz
    echo 'cat' ${f1} ${f2} ' > ' ${f3}
    cat ${f1} ${f2} > ${f3} &
    f1=${deep_data_path}${sample}/*${sample}*R2*fastq.gz
    f2=${shallow_data_path}${sample}/*${sample}*R2*fastq.gz
    f3=${target_data_path}fastq/${sample}/${sample}_R2.fastq.gz
    echo 'cat' ${f1} ${f2} ' > ' ${f3} 
    cat ${f1} ${f2} > ${f3} &
done
wait
echo merge done