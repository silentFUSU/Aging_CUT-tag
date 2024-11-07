raw_data=/mnt/transposon1/zhangyanxiaoLab/niuyuxiao/projects/liver_transposon/analysis/Dong_L_et_al_pri_hm_Ery_subclass/HiC/HiCPro_out/
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/Dong_L_et_al_pri_hm_Ery_subclass/
mkdir -p ${data_path}ValidPairs
pairix=/storage/zhangyanxiaoLab/suzhuojie/software/pairix/bin/
samples=(MergeAll_Ortho_HiC2_batch2_rep1 MergeAll_Ortho_HiC2_batch2_rep2 MergeAll_Pro_HiC2_batch2_rep1 MergeAll_Pro_HiC2_batch2_rep2)
chromsize=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/hg38/hg38.chrom.sizes
for sample in ${samples[@]}
do
    if ls ${raw_data}hic_results/data/${sample}/${sample}*.allValidPairs 1> /dev/null 2>&1; then  
        cp ${raw_data}hic_results/data/${sample}/${sample}*.allValidPairs ${data_path}ValidPairs/ &
    fi
    if ls ${raw_data}hic_results/data/${sample}*.allValidPairs.hic 1> /dev/null 2>&1; then  
        cp ${raw_data}hic_results/data/${sample}*.allValidPairs.hic ${data_path}juicer/ &
    fi
done
mkdir -p ${data_path}4DN_pairs
for sample in ${samples[@]}
do
    awk '{FS="\t";OFS="\t"} {print $1,$2,$3,$5,$6,$4,$7}' ${data_path}ValidPairs/${sample}.allValidPairs > ${data_path}4DN_pairs/${sample}.4DNpairs &
done

for sample in ${samples[@]}
do
    sort -k2,2 -k4,4 -k3,3n -k5,5n -S 10G ${data_path}4DN_pairs/${sample}.4DNpairs -T ${data_path}tmp/ > ${data_path}4DN_pairs/${sample}.4DNpairs.sort &
done 

for sample in ${samples[@]}
do
    rm ${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs
    ${pairix}bgzip ${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs.sort &
done

for sample in ${samples[@]}
do
    ${pairix}pairix -p pairs -f ${data_path}4DN_pairs/${sample}.4DNpairs.sort.gz &
done

for sample in ${samples[@]}
do
    python ~/software/pairsqc/pairsqc.py -M 8.4 -p ${data_path}4DN_pairs/${sample}.4DNpairs.sort.gz -c ${chromsize} -t P -O ${data_path}4DN_pairs/${sample} &
done

for sample in ${samples[@]}
do
    /usr/local/lib64/R/bin/Rscript ~/software/pairsqc/plot.r 4 ${data_path}4DN_pairs/${sample}_report/ &
done