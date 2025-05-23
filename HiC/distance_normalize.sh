tissue=$1
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
resolution=50000
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
mkdir -p ${data_path}${tissue}/ob_ex_matrix

if [[ "$tissue" == "ovary" || "$tissue" == "mammarygland" || "$tissue" == "uterus" ]]; then
    chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX)
else
    chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
fi

for sample in ${samples_for_tissue[@]}
do
    mkdir -p ${data_path}${tissue}/ob_ex_matrix/${sample}
    tag=${data_path}/${tissue}/compartment/homer_compartment/tagDir/${sample}/
    for chr in ${chromosomes[@]}
    do
        analyzeHiC ${tag} -chr ${chr} -res ${resolution} -distNorm -nolog > ${data_path}${tissue}/ob_ex_matrix/${sample}/${sample}_${resolution}_${chr}_ob_ex_Matrix.txt
    done
done