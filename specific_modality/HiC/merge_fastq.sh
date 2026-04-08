deep_data_path=/storage/zhangyanxiaoLab/fastq/2025/2025-06-30-Nuohe-DYQ/
shallow_data_path=/storage/zhangyanxiaoLab/fastq/2025/2025-06-18-Nuohe-DYQ/
target_data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data_transposon2/20250704_DYQ_HiC/
# samples=(WJH_stomach_96 WJH_stomach_99 WJH_stomach_108 WJH_stomach_112 WJH_thymus_230 WJH_thymus_240 WJH_thymus_241 WJH_thymus_245)
samples=(DYQ189 DYQ190)
mkdir -p ${target_data_path}
mkdir -p ${target_data_path}fastq/

for sample in ${samples[@]}
do
    mkdir -p ${target_data_path}fastq/
    mkdir -p ${target_data_path}fastq/${sample}
    f1=${deep_data_path}/*${sample}*/*${sample}*_1.fq.gz
    f2=${shallow_data_path}/*${sample}*/*${sample}*_1.fq.gz
    f3=${target_data_path}fastq/${sample}/${sample}_R1.fastq.gz
    echo 'cat' ${f1} ${f2} ' > ' ${f3}
    cat ${f1} ${f2} > ${f3} &
    f1=${deep_data_path}/*${sample}*/*${sample}*_2.fq.gz
    f2=${shallow_data_path}/*${sample}*/*${sample}*_2.fq.gz
    f3=${target_data_path}fastq/${sample}/${sample}_R2.fastq.gz
    echo 'cat' ${f1} ${f2} ' > ' ${f3} 
    cat ${f1} ${f2} > ${f3} &
done

# for sample in ${samples[@]}
# do
#     mkdir -p ${target_data_path}fastq/
#     mkdir -p ${target_data_path}fastq/${sample}
#     # f1=${deep_data_path}${sample}/*${sample}*R1*fastq.gz
#     # f2=${shallow_data_path}${sample}/*${sample}*R1*fastq.gz
#     f1=${deep_data_path}*${sample}*R1*fastq.gz
#     f2=${shallow_data_path}*${sample}*R1*fastq.gz
#     f3=${target_data_path}fastq/${sample}/${sample}_R1.fastq.gz
#     echo 'cat' ${f1} ${f2} ' > ' ${f3}
#     cat ${f1} ${f2} > ${f3} &
#     # f1=${deep_data_path}${sample}/*${sample}*R2*fastq.gz
#     # f2=${shallow_data_path}${sample}/*${sample}*R2*fastq.gz
#     f1=${deep_data_path}*${sample}*R2*fastq.gz
#     f2=${shallow_data_path}*${sample}*R2*fastq.gz
#     f3=${target_data_path}fastq/${sample}/${sample}_R2.fastq.gz
#     echo 'cat' ${f1} ${f2} ' > ' ${f3} 
#     cat ${f1} ${f2} > ${f3} &
# done
wait
echo merge done