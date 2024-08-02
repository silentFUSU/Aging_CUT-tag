data_path=$1
# data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/GSE119171_methyl_hic/fastq/
software_path=/storage/zhangyanxiaoLab/suzhuojie/software/
ref_data=/storage/zhangyanxiaoLab/suzhuojie/ref_data/
MARKDUP="/storage/zhangyanxiaoLab/share/Pipelines/atac-seq-pipeline-snakemake/dependencies/picard.jar MarkDuplicates"
sample=$2
ref=$3
# ref=mm10
# sample=SRR7770796
f1=${data_path}${sample}/*1.f*q.gz
f2=${data_path}${sample}/*2.f*q.gz
cd ${data_path}${sample}
mkdir ${data_path}${sample}/bam
mkdir ${data_path}${sample}/tmp/
# bismark_genome_preparation --path_to_aligner /usr/bin/ --verbose /storage/zhangyanxiaoLab/suzhuojie/ref_data/${ref}/
java -Xmx20G -Djava.library.path=${software_path}bisulfitehic/jbwa/src/main/native/ \
    -cp "${software_path}bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:${software_path}bisulfitehic/jbwa/jbwa.jar" main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem \
    ${ref_data}${ref}/${ref}.fa ${data_path}${sample}/bam/${sample}.bam ${f1} ${f2} -t 20 -rgId ${sample} -rgSm ${sample} -outputMateDiffChr -buffer 100000 

samtools sort --threads 10 -T ${sample}.name -n ${data_path}${sample}/bam/${sample}.bam | samtools fixmate -m --threads 10 - - | samtools sort --threads 10 -T ${data_path}${sample}/bam/${sample}.cor - | samtools markdup -T ${data_path}${sample}/bam/${sample}.mdups --threads 10 - - | samtools calmd --threads 10 -b - ${ref_data}${ref}/${ref}.fa 2>/dev/null > ${data_path}${sample}/bam/${sample}.calmd.bam
samtools index -@ 10 ${data_path}${sample}/bam/${sample}.calmd.bam

java -Xmx12G -jar ${MARKDUP} TMP_DIR=${data_path}${sample}/tmp/ INPUT=${data_path}${sample}/bam/${sample}.calmd.bam OUTPUT=${data_path}${sample}/bam/${sample}.calmd.nodup.bam METRICS_FILE=${data_path}${sample}/bam/${sample}.dup.qc VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true
samtools index -@ 10 ${data_path}${sample}/bam/${sample}.calmd.nodup.bam

mkdir ${data_path}${sample}/vcf
# java -Xmx5G -jar ${software_path}Bis-tools/Bis-SNP/BisSNP-1.0.0.jar -R /storage/zhangyanxiaoLab/suzhuojie/ref_data/${ref}/${ref}.fa \
#     -I  ${data_path}${sample}/bam/${sample}.calmd.bam -D /storage/zhangyanxiaoLab/suzhuojie/ref_data/dbsnp/dbsnp_138.b37.vcf.gz \
#     -T BisulfiteGenotyper -vfn1 ${data_path}${sample}/vcf/${sample}.calmd.cpg.default.raw.vcf -vfn2 ${data_path}${sample}/vcf/${sample}.calmd.snp.default.raw.vcf -C CG,1 -C CH,1 \
#     -out_modes DEFAULT_FOR_TCGA -stand_call_conf 20 -nt 30 -minConv 1 -vcfCache 1000000 -mmq 30 -mbq 5 


java -Xmx5G -jar ${software_path}Bis-tools/Bis-SNP/BisSNP-1.0.0.jar -R /storage/zhangyanxiaoLab/suzhuojie/ref_data/${ref}/${ref}.fa \
    -I  ${data_path}${sample}/bam/${sample}.calmd.nodup.bam -T BisulfiteGenotyper -vfn1 ${data_path}${sample}/vcf/${sample}.calmd.nodup.cpg.raw.vcf \
    -C CG,1 -C CH,1 -out_modes EMIT_ALL_CPG -stand_call_conf 20 -nt 30 -minConv 1 -vcfCache 1000000 -mmq 30 -mbq 5 

python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/methy_Hic/correct_vcf.py ${data_path}${sample}/vcf/${sample}.calmd.nodup.cpg.raw.vcf
perl ${software_path}Bis-tools/utils/sortByRefAndCor.pl --k 1 --c 2 --tmp ./ ${data_path}${sample}/vcf/corrected_${sample}.calmd.nodup.cpg.raw.vcf ${ref_data}${ref}/${ref}.fa.fai > ${data_path}${sample}/vcf/${sample}.calmd.nodup.cpg.raw.sort.vcf 
# perl ${software_path}Bis-tools/utils/sortByRefAndCor.pl --k 1 --c 2 --tmp ./ ${data_path}${sample}/${sample}.calmd.snp.raw.vcf ~/ref_data/${ref}/${ref}.fa.fai > ${data_path}${sample}/${sample}.calmd.snp.raw.sort.vcf 

# java -Xmx5G -jar ${software_path}Bis-tools/Bis-SNP/BisSNP-1.0.0.jar -R ${ref_data}${ref}/${ref}.fa -T VCFpostprocess -qual 20 -C CG -C CH -oldVcf ${data_path}${sample}/${sample}.calmd.cpg.raw.sort.vcf -snpVcf ${data_path}${sample}/${sample}.calmd.snp.raw.sort.vcf -newVcf ${data_path}${sample}/${sample}.calmd.cpg.filtered.sort.vcf -o ${data_path}${sample}/${sample}.calmd.cpg.filtered.sort.vcf.cpgSummary.txt -minCT 1
# java -Xmx5G -jar ${software_path}Bis-tools/Bis-SNP/BisSNP-1.0.0.jar -R ${ref_data}${ref}/${ref}.fa -T VCFpostprocess -qual 20 -C CG -C CH -oldVcf ${data_path}${sample}/${sample}.calmd.snp.raw.sort.vcf -snpVcf ${data_path}${sample}/${sample}.calmd.snp.raw.sort.vcf -newVcf ${data_path}${sample}/${sample}.calmd.snp.filtered.sort.vcf -o ${data_path}${sample}/${sample}.calmd.snp.filtered.sort.vcf.cpgSummary.txt 

# perl ${software_path}Bis-tools/utils/vcf2bedGraph.pl ${data_path}${sample}/${sample}.calmd.cpg.filtered.sort.vcf CG 
perl ${software_path}Bis-tools/utils/vcf2bed6plus2.pl --only_good_call ${data_path}${sample}/vcf/${sample}.calmd.nodup.cpg.raw.sort.vcf CG &
perl ${software_path}Bis-tools/utils/vcf2wig.pl ${data_path}${sample}/vcf/${sample}.calmd.nodup.cpg.raw.sort.vcf CG &
# perl ${software_path}Bis-tools/utils/vcf2coverage.pl  ${data_path}${sample}/${sample}.calmd.cpg.filtered.sort.vcf CG 
# ${software_path}Bis-tools/External_tools/ucsc_tools/bedGraphToBigWig ${data_path}${sample}/${sample}.calmd.cpg.filtered.sort.BisSNP-1.0.0.CG.coverage.bedgraph ${ref_data}/${ref}/${ref}.chrom.sizes ${data_path}${sample}/${sample}.calmd.cpg.filtered.sort.CG.bw 
mkdir ${data_path}${sample}/HiC
python ${software_path}bisulfitehic/src/python/sam2juicer_new.py -s  ${data_path}${sample}/bam/${sample}.calmd.nodup.bam -f ${ref_data}/${ref}/${ref}_DpnII.txt | gzip -fc > ${data_path}${sample}/HiC/${sample}.hic.nodup.txt.gz
zcat ${data_path}${sample}/HiC/${sample}.hic.nodup.txt.gz | /storage/zhangyanxiaoLab/share/Pipelines/hic-pipeline/lib/sort -k2,2 -k6,6 -k3,3n -k7,7n -u | gzip > ${data_path}${sample}/HiC/${sample}.hic.sortunique.nodup.txt.gz &

java -Xmx20G -jar ${software_path}juicer/CPU/common/juicer_tools.jar pre -d -q 30 -f ${ref_data}/${ref}/${ref}_DpnII.txt ${data_path}${sample}/HiC/${sample}.hic.nodup.txt.gz ${data_path}${sample}/HiC/${sample}.nodup.mapQ30.hic ${ref} &
# /usr/local/lib64/R/bin/Rscript ${software_path}bisulfitehic/src/R/call_topdom_on_hic.R ${data_path}${sample}/HiC/${sample}.nodup.mapQ30.hic ${data_path}${sample}/HiC/${sample}.nodup.mapQ30.topdom_default.tsv &
# python ${software_path}mustache/mustache/mustache.py -f ${data_path}${sample}/${sample}.mapQ30.hic -ch 22 -p 5 -r 25000 -cz ${ref_data}/${ref}/${ref}.chrom.sizes -norm KR -d 5000000 -o ${data_path}${sample}/${sample}.chr22.mustache.tsv

wait

fastq=${f1}
if [[ ${fastq: -3} == ".gz" ]]; then  
    Read="zcat"  
elif [[ ${fastq: -4} == ".bz2" ]]; then  
    Read="bzcat"  
else  
    Read="cat"  
fi  
fa=$($Read $fastq | wc -l)  
fa=$((fa/4))  
vpSortUnique=$(zcat  ${data_path}${sample}/HiC/${sample}.hic.sortunique.nodup.txt.gz | wc -l | cut -f 1 -d' ') 
# zcat ${data_path}${sample}/${sample}.hic.txt.gz | /storage/zhangyanxiaoLab/share/Pipelines/hic-pipeline/lib/sort -k2,2 -k6,6 -k3,3n -k7,7n -u | gzip > ${data_path}${sample}/${sample}.hic.nodup.txt.gz
vpNodup=$(zcat  ${data_path}${sample}/HiC/${sample}.hic.nodup.txt.gz | wc -l | cut -f 1 -d' ') 
cis_pair_15k=$(zcat ${data_path}${sample}/HiC/${sample}.hic.nodup.txt.gz | awk -v FS='\t' -v OFS='\t' '{if ($2 == $6) { m++; if ($7-$3 >= 15000 || $7-$3 <= -15000) { n++ } }} END { if (m == 0) { m = 0 } if (n == 0) { n = 0 } print m, n }')   
cis_pair_10k=$(zcat ${data_path}${sample}/HiC/${sample}.hic.nodup.txt.gz | awk -v FS='\t' -v OFS='\t' '{if ($2 == $6) { m++; if ($7-$3 >= 10000 || $7-$3 <= -10000) { n++ } }} END { if (m == 0) { m = 0 } if (n == 0) { n = 0 } print  n }')   
cis_pair_1k=$(zcat ${data_path}${sample}/HiC/${sample}.hic.nodup.txt.gz | awk -v FS='\t' -v OFS='\t' '{if ($2 == $6) { m++; if ($7-$3 >= 1000 || $7-$3 <= -1000) { n++ } }} END { if (m == 0) { m = 0 } if (n == 0) { n = 0 } print  n }')   
cis_pair_20k=$(zcat ${data_path}${sample}/HiC/${sample}.hic.nodup.txt.gz | awk -v FS='\t' -v OFS='\t' '{if ($2 == $6) { m++; if ($7-$3 >= 20000 || $7-$3 <= -20000) { n++ } }} END { if (m == 0) { m = 0 } if (n == 0) { n = 0 } print  n }')   
echo -e "Sample\tRaw\tVP_uniq\tVP_SortUnique\tVP_cis\tVP_cis15k\tVP_cis10k\tVP_cis1k\tVP_cis20k" > "${data_path}${sample}/HiC/${sample}.mapping_summary.txt"
echo -e "$sample\t$fa\t$vpNodup\t$vpSortUnique\t$cis_pair_15k\t$cis_pair_10k\t$cis_pair_1k\t$cis_pair_20k" >> "${data_path}${sample}/HiC/${sample}.mapping_summary.txt"

grep '^#' ${data_path}${sample}/vcf/${sample}.calmd.nodup.cpg.raw.sort.vcf  >  ${data_path}${sample}/vcf/tmp.${sample}.calmd.nodup.cpg.raw.sort.header.vcf 
grep -v '^#' ${data_path}${sample}/vcf/${sample}.calmd.nodup.cpg.raw.sort.vcf   | grep '^chrL' > ${data_path}${sample}/vcf/tmp.${sample}.chrL.calmd.nodup.cpg.raw.sort.vcf  
cat ${data_path}${sample}/vcf/tmp.${sample}.calmd.nodup.cpg.raw.sort.header.vcf  ${data_path}${sample}/vcf/tmp.${sample}.chrL.calmd.nodup.cpg.raw.sort.vcf \
  > ${data_path}${sample}/vcf/${sample}.chrL.calmd.nodup.cpg.raw.sort.vcf
rm ${data_path}${sample}/vcf/tmp*

/usr/local/lib64/R/bin/Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/methy_Hic/lambda_DNA_conversion.R ${data_path}${sample}/vcf/${sample}.chrL.calmd.nodup.cpg.raw.sort.vcf

mkdir ${data_path}${sample}/qc
mv ${data_path}${sample}/vcf/*MethySummarizeList.txt ${data_path}${sample}/qc
mv ${data_path}${sample}/HiC/*summary.txt ${data_path}${sample}/qc

# perl ~/software/Bis-tools/Bis-SNP/bissnp_easy_usage.pl --mmq 30 --nt 30 --mem 5 ~/software/Bis-tools/Bis-SNP/BisSNP-1.0.0.jar ${sample}.calmd.bam ${ref_data}${ref}/${ref}.fa  ${ref_data}dbsnp/dbsnp_138.b37.vcf.gz
# perl /your/path/to/Bis-SNP/bissnp_easy_usage.pl --use_bad_mates --mmq 30 --nt 30 --mem 5 /your/path/to/Bis-SNP/BisSNP-0.90.jar testOut.calmd.bam /your/path/to/reference_genome/human_g1k_v37.fa /your/path/to/dbsnp/dbsnp_138.b37.vcf
# /storage/zhangyanxiaoLab/suzhuojie/software/TrimGalore-0.6.10/trim_galore  ${f1} ${f2} --paired  \
#     --cores 4 --output_dir ${data_path}fastq_trim --basename {wildcards.sample} 2> {output.log}