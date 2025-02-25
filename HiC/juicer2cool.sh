source /storage/zhangyanxiaoLab/suzhuojie/miniconda3/etc/profile.d/conda.sh
conda activate hicexplorer
tissue=$1
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
resolution=50000
mkdir -p ${data_path}${tissue}/cool/
for sample in ${samples[@]}
do
    hic2cool convert ${data_path}${tissue}/juicer/${sample}.allValidPairs.hic ${data_path}${tissue}/cool/${sample}_${resolution}.cool -r ${resolution} &
done
wait 
echo all done