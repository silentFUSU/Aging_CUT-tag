tissue=ileum
raw_data=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240902_DYQ_WGBS/
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/WGBS_search_table.csv
samples=$(awk -F',' -v t="$tissue" -v a="$antibody" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  

for sample in ${samples[@]}
do
    mv ${raw_data}bed/${sample}* ${data_path}${tissue}/bdg/
    mv ${raw_data}bw/${sample}* ${data_path}${tissue}/bw/
done