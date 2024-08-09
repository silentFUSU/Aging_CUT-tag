cd /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240804_WGBS/compare/XX315/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/miniconda3/envs/snakemake/lib/
samtools index WGBS.bam
# Convert BAM to SAM  
samtools view -h juicer_meth.bam > juicer_meth.sam &  
samtools view -h WGBS.bam > WGBS.sam &
samtools view -h WGBS_local.bam  > WGBS_local.sam  &

# Extract read names  
awk '{if($1 !~ /^@/) print $1}' juicer_meth.sam > juicer_meth_readnames.txt  &
awk '{if($1 !~ /^@/) print $1}' WGBS.sam > WGBS_readnames.txt  &
# awk '{if($1 !~ /^@/) print $1}' WGBS_local.sam > WGBS_local_readnames.txt &

awk -F_ '{print $1}' WGBS_readnames.txt > WGBS_readnames_trimmed.txt  
# awk -F_ '{print $1}' WGBS_local_readnames.txt > WGBS_local_readnames_trimmed.txt  

sort -u juicer_meth_readnames.txt > juicer_meth_readnames_sort.txt  &
sort -u WGBS_readnames_trimmed.txt   > WGBS_readnames_sorted.txt  &
# sort WGBS_local_readnames_trimmed.txt >  WGBS_local_readnames_sorted.txt &  



# comm -23 juicer_meth_readnames_sort.txt WGBS_readnames_sorted.txt   > juicer_meth_readnames_unique.txt  

grep '^@' juicer_meth.sam > juicer_meth_header.sam 
grep -Ff juicer_meth_unique_readnames_sort.txt  juicer_meth.sam > juicer_meth_unique_reads_body.sam 
cat juicer_meth_header.sam  juicer_meth_unique_reads_body.sam  > juicer_meth_unique_reads.sam  


samtools view -bS juicer_meth_unique_reads.sam  | samtools sort -o juicer_meth_unique_reads.bam  
samtools index juicer_meth_unique_reads.bam  


/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant -q 30 -p 20 -d 1  \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.fa \
    juicer_meth_unique_reads.bam  

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
' juicer_meth_unique_reads_CpG.bedGraph

grep -Ff all_150M.txt juicer_meth.sam > juicer_unique_150M_body.sam 
cat juicer_meth_header.sam  juicer_unique_150M_body.sam    > 150M_test/juicer_unique_150M.sam 
samtools view -bS 150M_test/juicer_unique_150M.sam     | samtools sort -o 150M_test/juicer_unique_150M.bam 
samtools index 150M_test/juicer_unique_150M.bam 
/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant -q 30 -p 20 -d 1  \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.fa \
    150M_test/juicer_unique_150M.bam  
     
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
' 150M_test/juicer_unique_150M_CpG.bedGraph


grep -Ff juicer_unique_150M_short_distance_same_chr_dif_strand.txt juicer_meth.sam > juicer_unique_150M_short_distance_same_chr_dif_strand_body.sam 
cat juicer_meth_header.sam  juicer_unique_150M_short_distance_same_chr_dif_strand_body.sam   > 150M_test/juicer_unique_150M_short_distance_same_chr_dif_strand.sam 
samtools view -bS juicer_unique_150M_short_distance_same_chr_dif_strand.sam    | samtools sort -o juicer_unique_150M_short_distance_same_chr_dif_strand.bam 
samtools index juicer_unique_150M_short_distance_same_chr_dif_strand.bam   

/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant -q 30 -p 20 -d 1  \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.fa \
    juicer_unique_150M_short_distance_same_chr_dif_strand.bam   

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
' juicer_unique_150M_short_distance_same_chr_dif_strand_CpG.bedGraph #0.787001


grep -Ff juicer_unique_150M_dif_chr.txt juicer_meth.sam > juicer_unique_150M_dif_chr_body.sam 
cat juicer_meth_header.sam  juicer_unique_150M_dif_chr_body.sam    > 150M_test/juicer_unique_150M_dif_chr.sam 
samtools view -bS 150M_test/juicer_unique_150M_dif_chr.sam     | samtools sort -o 150M_test/juicer_unique_150M_dif_chr.bam 
samtools index 150M_test/juicer_unique_150M_dif_chr.bam 
/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant -q 30 -p 20 -d 1  \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.fa \
    150M_test/juicer_unique_150M_dif_chr.bam  
     
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
' 150M_test/juicer_unique_150M_dif_chr_CpG.bedGraph

grep -Ff juicer_unique_150M_same_chr_same_strand.txt juicer_meth.sam > juicer_unique_150M_same_chr_same_strand_body.sam 
cat juicer_meth_header.sam  juicer_unique_150M_same_chr_same_strand_body.sam    > 150M_test/juicer_unique_150M_same_chr_same_strand.sam 
samtools view -bS 150M_test/juicer_unique_150M_same_chr_same_strand.sam     | samtools sort -o 150M_test/juicer_unique_150M_same_chr_same_strand.bam 
samtools index 150M_test/juicer_unique_150M_same_chr_same_strand.bam 
/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant -q 30 -p 20 -d 1  \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.fa \
    150M_test/juicer_unique_150M_same_chr_same_strand.bam  
     
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
' 150M_test/juicer_unique_150M_same_chr_same_strand_CpG.bedGraph

grep -Ff juicer_unique_150M_same_chr_dif_strand_larger_500bp.txt juicer_meth.sam > juicer_unique_150M_same_chr_dif_strand_larger_500bp_body.sam 
cat juicer_meth_header.sam  juicer_unique_150M_same_chr_dif_strand_larger_500bp_body.sam    > 150M_test/juicer_unique_150M_same_chr_dif_strand_larger_500bp.sam 
samtools view -bS 150M_test/juicer_unique_150M_same_chr_dif_strand_larger_500bp.sam     | samtools sort -o 150M_test/juicer_unique_150M_same_chr_dif_strand_larger_500bp.bam 
samtools index 150M_test/juicer_unique_150M_same_chr_dif_strand_larger_500bp.bam 
/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant -q 30 -p 20 -d 1  \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.fa \
    150M_test/juicer_unique_150M_same_chr_dif_strand_larger_500bp.bam  
     
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
' 150M_test/juicer_unique_150M_same_chr_dif_strand_larger_500bp_CpG.bedGraph

grep -Ff juicer_unique_150M_same_chr_dif_strand_smaller_500bp.txt juicer_meth.sam > juicer_unique_150M_same_chr_dif_strand_smaller_500bp_body.sam 
cat juicer_meth_header.sam  juicer_unique_150M_same_chr_dif_strand_smaller_500bp_body.sam    > 150M_test/juicer_unique_150M_same_chr_dif_strand_smaller_500bp.sam 
samtools view -bS 150M_test/juicer_unique_150M_same_chr_dif_strand_smaller_500bp.sam     | samtools sort -o 150M_test/juicer_unique_150M_same_chr_dif_strand_smaller_500bp.bam 
samtools index 150M_test/juicer_unique_150M_same_chr_dif_strand_smaller_500bp.bam 
/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant -q 30 -p 20 -d 1  \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.fa \
    150M_test/juicer_unique_150M_same_chr_dif_strand_smaller_500bp.bam  
     
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
' 150M_test/juicer_unique_150M_same_chr_dif_strand_smaller_500bp_CpG.bedGraph


grep -Ff juicer_unique_Mpercent_smaller_90.txt juicer_meth.sam > juicer_unique_Mpercent_smaller_90_body.sam 
cat juicer_meth_header.sam  juicer_unique_Mpercent_smaller_90_body.sam    > juicer_unique_Mpercent_smaller_90.sam 
samtools view -bS juicer_unique_Mpercent_smaller_90.sam     | samtools sort -o juicer_unique_Mpercent_smaller_90.bam 
samtools index juicer_unique_Mpercent_smaller_90.bam 
/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant -q 30 -p 20 -d 1  \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.fa \
    juicer_unique_Mpercent_smaller_90.bam  
     
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
' juicer_unique_Mpercent_smaller_90_CpG.bedGraph

grep -Ff juicer_unique_Mpercent_larger90_smaller100.txt juicer_meth.sam > juicer_unique_Mpercent_larger90_smaller100_body.sam 
cat juicer_meth_header.sam  juicer_unique_Mpercent_larger90_smaller100_body.sam    > juicer_unique_Mpercent_larger90_smaller100.sam 
samtools view -bS juicer_unique_Mpercent_larger90_smaller100.sam     | samtools sort -o juicer_unique_Mpercent_larger90_smaller100.bam 
samtools index juicer_unique_Mpercent_larger90_smaller100.bam 
/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant -q 30 -p 20 -d 1  \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.fa \
    juicer_unique_Mpercent_larger90_smaller100.bam  
     
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
' juicer_unique_Mpercent_larger90_smaller100_CpG.bedGraph

grep -Ff juicer_unique_NOT_only_two_end.txt juicer_meth.sam > juicer_unique_NOT_only_two_end_body.sam 
cat juicer_meth_header.sam  juicer_unique_NOT_only_two_end_body.sam    > juicer_unique_NOT_only_two_end.sam 
samtools view -bS juicer_unique_NOT_only_two_end.sam     | samtools sort -o juicer_unique_NOT_only_two_end.bam 
samtools index juicer_unique_NOT_only_two_end.bam 
/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant -q 30 -p 20 -d 1  \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.fa \
    juicer_unique_NOT_only_two_end.bam  
     
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
' juicer_unique_NOT_only_two_end_CpG.bedGraph


cd /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240804_WGBS/compare/XX315_R1/
ln -s ../XX315/juicer_unique_150M_same_chr_same_strand.txt

samtools view -h XX315_R1_dedup_sorted.bam > XX315_R1_dedup_sorted.sam 
grep '^@' XX315_R1_dedup_sorted.sam  >  XX315_R1_dedup_sorted_header.sam 
grep -Ff juicer_unique_150M_same_chr_same_strand.txt  XX315_R1_dedup_sorted.sam  > XX315_R1_dedup_sorted_juicer_unique_150M_same_chr_same_strand_body.sam 
cat XX315_R1_dedup_sorted_header.sam   XX315_R1_dedup_sorted_juicer_unique_150M_same_chr_same_strand_body.sam  > XX315_R1_dedup_sorted_juicer_unique_150M_same_chr_same_strand.sam 
samtools view -bS XX315_R1_dedup_sorted_juicer_unique_150M_same_chr_same_strand.sam      | samtools sort -o XX315_R1_dedup_sorted_juicer_unique_150M_same_chr_same_strand.bam 
samtools index XX315_R1_dedup_sorted_juicer_unique_150M_same_chr_same_strand.bam #48912
/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant -q 30 -p 20 -d 1  \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.fa \
    XX315_R1_dedup_sorted_juicer_unique_150M_same_chr_same_strand.bam 

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
' XX315_R1_dedup_sorted_juicer_unique_150M_same_chr_same_strand_CpG.bedGraph #0.812585

cd /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240804_WGBS/compare/XX315_R2/
ln -s ../XX315/juicer_unique_150M_same_chr_same_strand.txt
ln -s /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240804_WGBS/WGBS_1M_single_end/XX315_R2/XX315_R2/bam/XX315_R2_dedup_sorted.bam
samtools view -h XX315_R2_dedup_sorted.bam > XX315_R2_dedup_sorted.sam 
grep '^@' XX315_R2_dedup_sorted.sam  >  XX315_R2_dedup_sorted_header.sam 
grep -Ff juicer_unique_150M_same_chr_same_strand.txt  XX315_R2_dedup_sorted.sam  > XX315_R2_dedup_sorted_juicer_unique_150M_same_chr_same_strand_body.sam 
cat XX315_R2_dedup_sorted_header.sam   XX315_R2_dedup_sorted_juicer_unique_150M_same_chr_same_strand_body.sam  > XX315_R2_dedup_sorted_juicer_unique_150M_same_chr_same_strand.sam 
samtools view -bS XX315_R2_dedup_sorted_juicer_unique_150M_same_chr_same_strand.sam      | samtools sort -o XX315_R2_dedup_sorted_juicer_unique_150M_same_chr_same_strand.bam 
samtools index XX315_R2_dedup_sorted_juicer_unique_150M_same_chr_same_strand.bam #48811
/storage/zhangyanxiaoLab/suzhuojie/software/MethylDackel/MethylDackel extract --keepSingleton --keepDiscordant -q 30 -p 20 -d 1  \
    /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.fa \
    XX315_R2_dedup_sorted_juicer_unique_150M_same_chr_same_strand.bam 

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
' XX315_R2_dedup_sorted_juicer_unique_150M_same_chr_same_strand_CpG.bedGraph #0.0206547

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
' merged_dedup_sort_CpG.bedGraph