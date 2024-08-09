data_path=/storage/zhangyanxiaoLab/zhangyanxiao/software/obsutil_linux_amd64_5.2.10/tmp/XX/20240804-1/
target_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240804_WGBS/bhmem_1M/
samples=(XX315 XX316)
for sample in ${samples[@]}
do
    mkdir ${target_path}${sample}
    mkdir ${target_path}${sample}/fastq
    zcat ${data_path}${sample}/*R1*.fastq.gz | head -n 4000000 | gzip > ${target_path}${sample}/fastq/${sample}_R1_001.fastq.gz &
    zcat ${data_path}${sample}/*R2*.fastq.gz | head -n 4000000 | gzip > ${target_path}${sample}/fastq/${sample}_R2_001.fastq.gz &
done

data_path=/mnt/transposon1/zhangyanxiaoLab/suzhuojie/project/Aging_CUT_TAG/public_data/GSE119171_NMethod_methylHiC/fastq/
target_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/GSE119171_methyl_hic/SRR7770813_WGBS/single_end/fastq/
samples=(SRR7770813)
for sample in ${samples[@]}
do
    zcat ${data_path}${sample}_1.fastq.gz | head -n 4000000 | gzip > ${target_path}${sample}_R1.fastq.gz &
    zcat ${data_path}${sample}_2.fastq.gz | head -n 4000000 | gzip > ${target_path}${sample}_R2.fastq.gz &
done

data_path=/storage/zhangyanxiaoLab/fastq/2024/2024-07-25-Jiangbei-YuLab/WJH_Mousebrain_C1_BSseq/
target_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/light_test/WJH_Mousebrain_C1_BSseq/WJH_Mousebrain_C1_BSseq_101_single_end-end-to-end/fastq/
samples=(WJH_Mousebrain_101_C1_Bsseq)
for sample in ${samples[@]}
do
    zcat ${data_path}${sample}*R1*.fastq.gz | head -n 4000000 | gzip > ${target_path}${sample}_R1.fastq.gz &
    zcat ${data_path}${sample}*R2*.fastq.gz | head -n 4000000 | gzip > ${target_path}${sample}_R2.fastq.gz &
done

data_path=/storage/zhangyanxiaoLab/zhangyanxiao/software/obsutil_linux_amd64_5.2.10/tmp/XX/20240804-1/XX315/
target_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240804_WGBS/WGBS_1M_single_end_trim_70bp/fastq/
sample=XX315
zcat ${data_path}${sample}*R2*.fastq.gz | head -n 4000000 | gzip > ${target_path}${sample}_R2.fastq.gz &
trim_galore --length 20 --stringency 3 --cores 8 --gzip \
    --clip_R1 70 \
    ${target_path}${sample}_R2.fastq.gz 2>/dev/null  
mv XX315_R2_trimmed.fq.gz XX315_R2_5_trimmed_70.fastq.gz
mv XX315_R2.fastq.gz_trimming_report.txt XX315_R2.fastq.gz_5_trimming_70_report.txt

trim_galore --length 20 --stringency 3 --cores 8 --gzip \
    --three_prime_clip_R1 70 \
    ${target_path}${sample}_R2.fastq.gz 2>/dev/null  

mv XX315_R2_trimmed.fq.gz XX315_R2_3_trimmed_70.fastq.gz
mv XX315_R2.fastq.gz_trimming_report.txt XX315_R2.fastq.gz_3_trimming_70_report.txt

data_path=/storage/zhangyanxiaoLab/zhangyanxiao/software/obsutil_linux_amd64_5.2.10/tmp/XX/20240804-1/XX315/
target_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240804_WGBS/WGBS_1M_single_end/XX315_R2_pbat/fastq/
sample=XX315
zcat ${data_path}${sample}*R2*.fastq.gz | head -n 4000000 | gzip > ${target_path}${sample}_R2.fastq.gz &


data_path=/storage/zhangyanxiaoLab/zhangyanxiao/software/obsutil_linux_amd64_5.2.10/tmp/XX/20240804-1/XX315/
target_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240804_WGBS/WGBS_1M_single_end_trim/fastq/
sample=XX315
zcat ${data_path}${sample}*R2*.fastq.gz | head -n 4000000 | gzip > ${target_path}${sample}_R2.fastq.gz &
trim_galore --length 20 --stringency 3 --cores 8 --gzip \
    --clip_R1 10 \
    ${target_path}${sample}_R2.fastq.gz 2>/dev/null  
mv XX315_R2_trimmed.fq.gz XX315_R2_5_trimmed_10.fastq.gz
mv XX315_R2.fastq.gz_trimming_report.txt XX315_R2.fastq.gz_5_trimming_10_report.txt