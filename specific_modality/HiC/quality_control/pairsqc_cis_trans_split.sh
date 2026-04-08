data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=lung
chromsize=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.normal.chrom.sizes
pairix=/storage/zhangyanxiaoLab/suzhuojie/software/pairix/bin/
mkdir -p ${data_path}${tissue}/4DN_pairs
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table") 

for sample in ${samples_for_tissue[@]}
do
    zcat "${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs.sort.gz" | awk -v cis="${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs.cis.sort" -v trans="${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs.trans.sort" '  
    {  
        if ($2 == $4) {  
            print $0 > cis  
        } else {  
            print $0 > trans  
        }  
    }'
    ${pairix}bgzip ${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs.cis.sort &
    ${pairix}bgzip ${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs.trans.sort &
    wait
    ${pairix}pairix -p pairs -f ${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs.cis.sort.gz &
done

for sample in ${samples_for_tissue[@]}
do
    python ~/software/pairsqc/pairsqc.py -M 8.3 -p ${data_path}${tissue}/4DN_pairs/${sample}.4DNpairs.cis.sort.gz -c ${chromsize} -t P -O ${data_path}${tissue}/4DN_pairs/${sample}_cis &
done

for sample in ${samples_for_tissue[@]}
do
    /usr/local/lib64/R/bin/Rscript ~/software/pairsqc/plot.r 4 ${data_path}${tissue}/4DN_pairs/${sample}_cis_report/ &
done