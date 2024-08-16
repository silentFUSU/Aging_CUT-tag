data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240804_WGBS/biscuit_WGBS_1M/
samples=(XX316)
# biscuit index ~/software/juicer/references/mm10_lambda/mm10_lambda.fa
for sample in ${samples[@]}
do
    mkdir ${data_path}${sample}/bam
    libraryName=WGBS_library
    rg="@RG\\tID:${sample}\\tSM:${sample}\\tPL:LS454\\tLB:${libraryName}"
    biscuit align -@ 16 -R $rg ~/software/juicer/references/mm10_lambda/mm10_lambda.fa ${data_path}${sample}/fastq/*R1*gz ${data_path}${sample}/fastq/*R2*gz | \
        dupsifter  ~/software/juicer/references/mm10_lambda/mm10_lambda.fa  | \
        samtools sort -@ 8 -o ${data_path}${sample}/bam/${sample}.bam -O BAM - 
    samtools index ${data_path}${sample}/bam/${sample}.bam 
    mkdir ${data_path}${sample}/${sample}_report
    cd ${data_path}${sample}/${sample}_report
    biscuit qc  ~/software/juicer/references/mm10_lambda/mm10_lambda.fa ${data_path}${sample}/bam/${sample}.bam ${sample}
    multiqc .
done

samtools view -b -q 40 ${data_path}${sample}/bam/${sample}.bam > ${data_path}${sample}/bam/${sample}_filter.bam
samtools sort  ${data_path}${sample}/bam/${sample}_filter.bam -o  ${data_path}${sample}/bam/${sample}_filter_sort.bam
mkdir ${data_path}${sample}/${sample}_filter_sort_report
cd ${data_path}${sample}/${sample}_filter_sort_report
biscuit qc  ~/software/juicer/references/mm10_lambda/mm10_lambda.fa ${data_path}${sample}/bam/${sample}_filter_sort.bam ${sample}

samples=(XX315 XX316)
for sample in ${samples[@]}
do
    mkdir ${data_path}${sample}/bam
    libraryName=WGBS_library
    rg="@RG\\tID:${sample}\\tSM:${sample}\\tPL:LS454\\tLB:${libraryName}"
    cd ${data_path}${sample}/fastq
    trim_galore --length 20 --stringency 3 --cores 8 --gzip \
        --clip_R1 10 \
        ${data_path}${sample}/fastq/${sample}_R2.fastq.gz 2>/dev/null
    biscuit align -@ 16 -R $rg ~/software/juicer/references/mm10_lambda/mm10_lambda.fa ${data_path}${sample}/fastq/${sample}_R2_trimmed.fq.gz | \
        dupsifter  ~/software/juicer/references/mm10_lambda/mm10_lambda.fa  | \
        samtools sort -@ 8 -o ${data_path}${sample}/bam/${sample}_R2.bam -O BAM - 
    samtools index ${data_path}${sample}/bam/${sample}_R2.bam 
done