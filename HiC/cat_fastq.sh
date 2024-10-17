raw_data=/storage/zhangyanxiaoLab/zhangyanxiao/software/obsutil_linux_amd64_5.2.10/tmp/YuLab/20241012/
target_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20241012_WJH_nomerge_HiC/fastq/
samples=(WJH-106-Lung)
for sample in ${samples[@]}
do
    echo $sample
    mkdir ${target_path}${sample}
    ln -s ${raw_data}${sample}*/*${sample}*_1*.gz ${target_path}${sample}/${sample}_R1.fastq.gz
    ln -s ${raw_data}${sample}*/*${sample}*_2*.gz ${target_path}${sample}/${sample}_R2.fastq.gz
done

raw_data=/storage/zhangyanxiaoLab/zhangyanxiao/software/obsutil_linux_amd64_5.2.10/tmp/YuLab/20240925-1/WJH-Ileum-test2/
target_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240925_WJH_HiC/fastq/
samples=(WJH-103-Ileum WJH-106-Ileum WJH-109-Ileum)
for sample in ${samples[@]}
do
    echo $sample
    mkdir ${target_path}${sample}
    ln -s ${raw_data}${sample}*R1*.gz ${target_path}${sample}/${sample}_R1.fastq.gz
    ln -s ${raw_data}${sample}*R2*.gz ${target_path}${sample}/${sample}_R2.fastq.gz
done

raw_data=/storage/zhangyanxiaoLab/zhangyanxiao/software/obsutil_linux_amd64_5.2.10/tmp/YuLab/20240918/
target_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240918_WJH_pipeline_test_HiC/fastq/
samples=(WJH-100-Lung WJH-103-Lung WJH-106-Lung WJH-109-Lung)
for sample in ${samples[@]}
do
    echo $sample
    mkdir ${target_path}${sample}
    zcat ${raw_data}${sample}*/${sample}*R1*.gz | head -n 4000000 | gzip > ${target_path}${sample}/${sample}_R1.fastq.gz &
    zcat ${raw_data}${sample}*/${sample}*R2*.gz | head -n 4000000 | gzip > ${target_path}${sample}/${sample}_R2.fastq.gz &
done

for sample in ${samples[@]}
do
    mkdir ${target_path}${sample}
    mkdir ${target_path}${sample}/fastq
    ln -s ${raw_data}${sample}*/${sample}*R1*.gz ${target_path}${sample}/fastq/${sample}_R1.fastq.gz
    ln -s ${raw_data}${sample}*/${sample}*R2*.gz ${target_path}${sample}/fastq/${sample}_R2.fastq.gz
done