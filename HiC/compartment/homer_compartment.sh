data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=$1
res=50000
mkdir -p ${data_path}${tissue}/homer_compartment
mkdir -p ${data_path}${tissue}/homer_compartment/allValidPair_homer
mkdir -p ${data_path}${tissue}/homer_compartment/tagDir
mkdir -p ${data_path}${tissue}/homer_compartment/PC1
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  

# samples_for_tissue=(WJH-100-Lung WJH-103-Lung WJH-106-Lung WJH-109-Lung)
for sample in ${samples_for_tissue[@]}
do
    valid_pairs=${data_path}${tissue}/ValidPairs/${sample}.allValidPairs
    cut -d $'\t' -f 1-7 ${valid_pairs} > ${data_path}${tissue}/homer_compartment/allValidPair_homer/${sample}.allValidPairs.homer
    makeTagDirectory ${data_path}${tissue}/homer_compartment/tagDir/${sample} -format HiCsummary \
        ${data_path}${tissue}/homer_compartment/allValidPair_homer/${sample}.allValidPairs.homer
    
    runHiCpca.pl ${data_path}${tissue}/homer_compartment/PC1/${sample}_${res} ${data_path}${tissue}/homer_compartment/tagDir/${sample} -res ${res} -cpu 1 -genome mm10
done

# conditions=(young old)
# for condition in ${conditions[@]}
# do
#     valid_pairs=${data_path}${tissue}/ValidPairs/${condition}_combined.allValidPairs
#     cut -d $'\t' -f 1-7 ${valid_pairs} > ${data_path}${tissue}/homer_compartment/allValidPair_homer/${condition}_combined.allValidPairs.homer
#     makeTagDirectory ${data_path}${tissue}/homer_compartment/tagDir/${condition}_combined -format HiCsummary \
#         ${data_path}${tissue}/homer_compartment/allValidPair_homer/${condition}_combined.allValidPairs.homer
    
#     runHiCpca.pl ${data_path}${tissue}/homer_compartment/PC1/${condition}_combined ${data_path}${tissue}/homer_compartment/tagDir/${condition}_combined -res 20000 -cpu 10 -genome mm10
# done