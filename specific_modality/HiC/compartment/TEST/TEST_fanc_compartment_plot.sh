data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/HiC/
tissue=lung
resolution=1mb
# position=chr19:11000001-60000000
position=chr19
positions=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
mkdir -p ${result_path}${tissue}/
mkdir -p ${result_path}${tissue}/compartment/
for position in ${positions[@]}
do
    fancplot -o ${result_path}${tissue}/compartment/fanc_${tissue}_young_${resolution}_${position}.ab_and_ev.png ${position} \
        -p square ${data_path}${tissue}/compartment/fanc_young_${resolution}.ab \
        -c RdBu_r -vmin -1 -vmax 1 \
        -p line ${data_path}${tissue}/compartment/fanc_young_${resolution}.ev.txt

    fancplot -o ${result_path}${tissue}/compartment/fanc_${tissue}_old_${resolution}_${position}.ab_and_ev.png ${position} \
        -p square ${data_path}${tissue}/compartment/fanc_old_${resolution}.ab \
        -c RdBu_r -vmin -1 -vmax 1 \
        -p line ${data_path}${tissue}/compartment/fanc_old_${resolution}.ev.txt
done

search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
ref=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.fa
for sample in ${samples_for_tissue[@]}
do
    fancplot -o ${result_path}${tissue}/compartment/fanc_${tissue}_${sample}_${resolution}_${position}.ab_and_ev.png ${position} \
    -p square ${data_path}${tissue}/compartment/fanc_${sample}_${resolution}.ab \
    -c RdBu_r -vmin -1 -vmax 1 \
    -p line ${data_path}${tissue}/compartment/fanc_${sample}_${resolution}.ev.txt
done