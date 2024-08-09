/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant -q 30 -p 20 -d 1  \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.fa \
    /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/bhmem/WJH_Mousebrain_C1_BSseq_101/bam/WJH_Mousebrain_C1_BSseq_101.calmd.nodup.bam

/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant  --cytosine_report  --CHH --CHG -q 30 -p 20 -d 1 \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.fa \
    /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/bhmem/WJH_Mousebrain_C1_BSseq_101/bam/WJH_Mousebrain_C1_BSseq_101.calmd.nodup.bam

software_path=/storage/zhangyanxiaoLab/suzhuojie/software/
ref_data=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_methylHiC/
sample=WJH_Mousebrain_C1_BSseq_101
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/bhmem/
ref=mm10
java -Xmx5G -jar ${software_path}Bis-tools/Bis-SNP/BisSNP-1.0.0.jar -R ${ref_data}${ref}/${ref}.fa \
    -I  ${data_path}${sample}/bam/${sample}.calmd.nodup.bam -T BisulfiteGenotyper -vfn1 ${data_path}${sample}/vcf_minconv0/${sample}.calmd.nodup.cpg.raw.vcf \
    -C CG,1 -C CH,1 -out_modes EMIT_ALL_CPG -stand_call_conf 20 -nt 30 -minConv 0 -vcfCache 1000000 -mmq 30 -mbq 5 
python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/methy_Hic/correct_vcf.py ${data_path}${sample}/vcf_minconv0/${sample}.calmd.nodup.cpg.raw.vcf
perl ${software_path}Bis-tools/utils/sortByRefAndCor.pl --k 1 --c 2 --tmp ./ ${data_path}${sample}/vcf_minconv0/corrected_${sample}.calmd.nodup.cpg.raw.vcf ${ref_data}${ref}/${ref}.fa.fai > ${data_path}${sample}/vcf_minconv0/${sample}.calmd.nodup.cpg.raw.sort.vcf 
# perl ${software_path}Bis-tools/utils/vcf2bed6plus2.pl --only_good_call ${data_path}${sample}/vcf_minconv0/${sample}.calmd.nodup.cpg.raw.sort.vcf CG &
perl ${software_path}Bis-tools/utils/vcf2bed6plus2.pl --qual 0  --maxCov 50000 ${data_path}${sample}/vcf_minconv0/${sample}.calmd.nodup.cpg.raw.sort.vcf CG &

outputdir=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/lambda_test_all_data/WJH_Mousebrain_C1_BSseq/WJH_Mousebrain_C1_BSseq_101/aligned/
samtools sort -@ 8 ${outputdir}/merged_dedup.bam > ${outputdir}/merged_dedup_sort.bam
samtools index  ${outputdir}/merged_dedup_sort.bam
/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant  --cytosine_report  --CHH --CHG \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10_lambda/mm10_lambda.fa \
    ${outputdir}/merged_dedup_sort.bam

awk '  
BEGIN {  
    meth = 0;  
    unmeth = 0;  
}  
{  
    if ($6 == "CG" && ($4 + $5) > 0) {  
        if ($4 > $5) {  
            meth++;  
        } else if ($4 < $5) {  
            unmeth++;  
        }  
    }  
}  
END {  
    print "meth: " meth;  
    print "unmeth: " unmeth;  
}  
' ${outputdir}merged_dedup_sort.cytosine_report.txt

outputdir=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/lambda_test_all_data/WJH_Mousebrain_C1_BSseq/WJH_Mousebrain_C1_BSseq_101/aligned/
awk '  
{  
    column5_sum += $5;  
    column6_sum += $6;  
}  
END {  
    total_sum = column5_sum + column6_sum;  
    if (total_sum != 0) {  
        ratio = column5_sum / total_sum;  
        print "第五列的总和除以（第五列的总和+第六列的总和）为：" ratio;  
    } else {  
        print "第五列和第六列的总和都为0，无法计算比例。";  
    }  
}  
' ${outputdir}merged_dedup_sort_CpG.bedGraph