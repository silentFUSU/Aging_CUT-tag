cat /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/mm10.fa /storage/zhangyanxiaoLab/suzhuojie/ref_data/lambda/lambda.fa > /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/mm10_lambda.fa
rm /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/mm10.fa
mv /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/mm10_lambda.fa /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/mm10.fa
samtools faidx /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/mm10.fa
bismark_genome_preparation --path_to_aligner /usr/bin/ --verbose /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/ &
cut -f1,2 /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/mm10.fa.fai > /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/mm10.chrom.sizes
cd /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/
python /storage/zhangyanxiaoLab/suzhuojie/software/juicer/misc/generate_site_positions.py DpnII mm10 /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/mm10.fa
java -jar /storage/zhangyanxiaoLab/share/Pipelines/atac-seq-pipeline-snakemake/dependencies/picard.jar CreateSequenceDictionary REFERENCE=/storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/mm10.fa  OUTPUT=/storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/mm10.dict
wait
samtools faidx /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa
samtools faidx /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa
bwa index /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa
bwa index /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa