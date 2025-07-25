data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/HBEC/sample/DYQ181/
chromsize=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/hg38/hg38.normal.chrom.sizes
pairix=/storage/zhangyanxiaoLab/suzhuojie/software/pairix/bin/
mkdir -p ${data_path}4DN_pairs
mkdir -p ${data_path}tmp

samples=(DYQ181-1 DYQ181-2 DYQ181-3 DYQ181-4)

for sample in ${samples[@]}
do
    awk '{FS="\t";OFS="\t"} {print $1,$2,$3,$5,$6,$4,$7}' ${data_path}ValidPairs/${sample}.allValidPairs > ${data_path}4DN_pairs/${sample}.4DNpairs &
done
wait

for sample in ${samples[@]}
do
    sort -k2,2 -k4,4 -k3,3n -k5,5n -S 10G ${data_path}4DN_pairs/${sample}.4DNpairs -T ${data_path}tmp/ > ${data_path}4DN_pairs/${sample}.4DNpairs.sort &
done 
wait

for sample in ${samples[@]}
do
    rm ${data_path}4DN_pairs/${sample}.4DNpairs
    ${pairix}bgzip ${data_path}4DN_pairs/${sample}.4DNpairs.sort &
done
wait

for sample in ${samples[@]}
do
    ${pairix}pairix -p pairs -f ${data_path}4DN_pairs/${sample}.4DNpairs.sort.gz &
done
wait

for sample in ${samples[@]}
do
    python ~/software/pairsqc/pairsqc.py -M 8.3 -p ${data_path}4DN_pairs/${sample}.4DNpairs.sort.gz -c ${chromsize} -t P -O ${data_path}4DN_pairs/${sample} &
done
wait

for sample in ${samples[@]}
do
    /usr/local/lib64/R/bin/Rscript ~/software/pairsqc/plot.r 4 ${data_path}4DN_pairs/${sample}_report/ &
done
echo all done
# for line in open("/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/hg38/hg38.chrom.sizes"):  
#     if '\t' not in line:  
#         print(f"Warning: Incorrect format in line: {line}")  
#     else:  
#         chrom, size = line.strip().split('\t')  