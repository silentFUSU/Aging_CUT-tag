data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=lung
resolution=10000
juiceDir=/storage/zhangyanxiaoLab/suzhuojie/software/juicer/
mkdir -p ${data_path}lung/loop/
mkdir -p ${data_path}lung/loop/HiCCUPS
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
method=SCALE
chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)

for sample in ${samples_for_tissue[@]}
do
    mkdir -p ${data_path}${tissue}/loop/HiCCUPS/${sample}
    juicer_file=${data_path}${tissue}/juicer/${sample}.allValidPairs.hic
    for chr in ${chromosomes[@]}
    do
        mkdir -p ${data_path}${tissue}/loop/HiCCUPS/${sample}/${chr}
        ${juiceDir}/scripts/common/juicer_tools hiccups --cpu -r ${resolution} -k ${method} -f 0.1 -p 2 -i 8 -t 0.02,1.5,1.75,2 -d 40000 -c ${chr} ${juicer_file} ${data_path}${tissue}/loop/HiCCUPS/${sample}/${chr}
    done
done