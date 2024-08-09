data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240804_WGBS/biscuit_WGBS_1M/
samples=(XX315 XX316)
biscuit index ~/software/juicer/references/mm10_lambda/mm10_lambda.fa
for sample in ${samples[@]}
do
    mkdir ${data_path}${sample}/bam
    libraryName=WGBS_library

    rg="@RG\\tID:${sample}\\tSM:${sample}\\tPL:LS454\\tLB:${libraryName}"
    biscuit align -@ 16 -R $rg ~/software/juicer/references/mm10_lambda/mm10_lambda.fa ${data_path}${sample}/fastq/*R1*gz ${data_path}${sample}/fastq/*R2*gz | \
    dupsifter  ~/software/juicer/references/mm10_lambda/mm10_lambda.fa  | \
    samtools sort -@ 8 -o ${data_path}${sample}/bam/${sample}.bam -O BAM -
    samtools index ${data_path}${sample}/bam/${sample}.bam 
done