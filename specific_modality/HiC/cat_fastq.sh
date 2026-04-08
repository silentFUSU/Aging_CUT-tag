raw_data=/storage/zhangyanxiaoLab/suzhuojie/software/obsutil_linux_amd64_5.2.10/tmp/SZJ/20260331/
target_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data_transposon2/20260331_microC_DYQ/
# samples=(Cecum-103 Cecum-106 Cecum-109 WJH-Ileum-100 WJH-Ileum-103 WJH-Ileum-106 WJH-ileum-112 WJH-pancreas-107 WJH-pancreas-108 WJH-pancreas-98  WJH-thymus-230 WJH-thymus-240 WJH-thymus-241 WJH-thymus-245)
samples=(DYQ301)
mkdir -p ${target_path}
mkdir -p ${target_path}fastq/
for sample in ${samples[@]}
do
    echo $sample
    mkdir ${target_path}fastq/${sample}
    # ln -s ${raw_data}${sample}*/*${sample}*_R1*.gz ${target_path}fastq/${sample}/${sample}_R1.fastq.gz
    # ln -s ${raw_data}${sample}*/*${sample}*_R2*.gz ${target_path}fastq/${sample}/${sample}_R2.fastq.gz
    # ln -s ${raw_data}*${sample}*R1*.gz ${target_path}fastq/${sample}/${sample}_R1.fastq.gz
    # ln -s ${raw_data}*${sample}*R2*.gz ${target_path}fastq/${sample}/${sample}_R2.fastq.gz
    # mkdir ${target_path}${sample}
    # mkdir ${target_path}${sample}/fastq
    ln -s ${raw_data}${sample}*/*${sample}*_R1*.gz ${target_path}fastq/${sample}/${sample}_R1.fastq.gz
    ln -s ${raw_data}${sample}*/*${sample}*_R2*.gz ${target_path}fastq/${sample}/${sample}_R2.fastq.gz
done

raw_data=/storage/zhangyanxiaoLab/zhangyanxiao/software/obsutil_linux_amd64_5.2.10/tmp/YuLab/20240925-1/WJH-Ileum-test2/
target_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240925_WJH_HiC/fastq/
samples=(WJH-103-Ileum WJH-106-Ileum WJH-109-Ileum)
for sample in ${samples[@]}
do
    echo $sample
    mkdir ${target_path}${sample}
    ln -s ${raw_data}${sample}*_1*.gz ${target_path}${sample}/${sample}_R1.fastq.gz
    ln -s ${raw_data}${sample}*_2*.gz ${target_path}${sample}/${sample}_R2.fastq.gz
done

raw_data=/storage/zhangyanxiaoLab/fastq/2024/2024-11-11-Jiangbei-YuLab/20241127/
target_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data_transposon2/CB_HiC_test_nomerge/
samples=(WJH-Cerebellum-99 WJH-Cerebellum-108)
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
    ln -s ${raw_data}${sample}*/*${sample}*_1*.gz ${target_path}${sample}/fastq/${sample}_R1.fastq.gz
    ln -s ${raw_data}${sample}*/*${sample}*_2*.gz ${target_path}${sample}/fastq/${sample}_R2.fastq.gz
done