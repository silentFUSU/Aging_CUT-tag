data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=lung
resolution=10000
juiceDir=/storage/zhangyanxiaoLab/suzhuojie/software/juicer/
mkdir -p ${data_path}lung/TAD/arrowhead
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
method=SCALE

for sample in ${samples_for_tissue[@]}
do
    ${juiceDir}/scripts/common/juicer_tools arrowhead -k SCALE -r ${resolution} --threads 5 ${data_path}${tissue}/juicer/${sample}.allValidPairs.hic ${data_path}lung/TAD/arrowhead/${sample}_${resolution}_${method}_output
done