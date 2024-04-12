data_path=~/projects/Aging_CUT_Tag/data/raw_data/methy_HiC/20240411_hundep2_methylHic_data/
software_path=~/software/
ref_data=~/ref_data/
sample=WJH_MethylHic_10

f1=${data_path}${sample}/hundep2_MethylHic_10_S1_L002_R1_001.fastq.gz
f2=${data_path}${sample}/hundep2_MethylHic_10_S1_L002_R2_001.fastq.gz
# bismark_genome_preparation --path_to_aligner /usr/bin/ --verbose /storage/zhangyanxiaoLab/suzhuojie/ref_data/hg38/
nohup java -Xmx20G -Djava.library.path=/storage/zhangyanxiaoLab/suzhuojie/software/bisulfitehic/jbwa/src/main/native/ \
    -cp "/storage/zhangyanxiaoLab/suzhuojie/software/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/storage/zhangyanxiaoLab/suzhuojie/software/bisulfitehic/jbwa/jbwa.jar" main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem \
    ${ref_data}hg38/hg38.fa ${data_path}${sample}/${sample}.bam ${f1} ${f2} -t 5 -rgId ${sample} -rgSm ${sample} -outputMateDiffChr -buffer 100000 &


/storage/zhangyanxiaoLab/suzhuojie/software/TrimGalore-0.6.10/trim_galore  ${f1} ${f2} --paired  \
    --cores 4 --output_dir ${data_path}fastq_trim --basename {wildcards.sample} 2> {output.log}