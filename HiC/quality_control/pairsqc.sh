data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=lung
chromsize=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.normal.chrom.sizes
pairix=/storage/zhangyanxiaoLab/suzhuojie/software/pairix/bin/
mkdir -p ${data_path}${tissue}/4DN_pairs
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  

for sample in ${samples_for_tissue[@]}
do
    awk '{FS="\t";OFS="\t"} {print $1,$2,$3,$5,$6,$4,$7}' ${data_path}${tissue}/ValidPairs/${sample}.allValidPairs > ${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs &
done
wait

for sample in ${samples_for_tissue[@]}
do
    sort -k2,2 -k4,4 -k3,3n -k5,5n -S 10G ${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs -T ${data_path}${tissue}/tmp/ > ${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs.sort &
done 
wait

for sample in ${samples_for_tissue[@]}
do
    rm ${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs
    ${pairix}bgzip ${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs.sort &
done
wait

for sample in ${samples_for_tissue[@]}
do
    ${pairix}pairix -p pairs -f ${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs.sort.gz &
done
wait

for sample in ${samples_for_tissue[@]}
do
    python ~/software/pairsqc/pairsqc.py -M 8.3 -p ${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs.sort.gz -c ${chromsize} -t P -O ${data_path}${tissue}/4DN_pairs/${sample} &
done
wait

for sample in ${samples_for_tissue[@]}
do
    /usr/local/lib64/R/bin/Rscript ~/software/pairsqc/plot.r 4 ${data_path}${tissue}/4DN_pairs/${sample}_report/ &
done
echo all done
# for line in open("/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.chrom.sizes"):  
#     if '\t' not in line:  
#         print(f"Warning: Incorrect format in line: {line}")  
#     else:  
#         chrom, size = line.strip().split('\t')  