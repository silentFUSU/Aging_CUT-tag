data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=lung
resolution=10kb
mustacheDir=/storage/zhangyanxiaoLab/suzhuojie/software/mustache/mustache/
mkdir -p ${data_path}lung/loop/
mkdir -p ${data_path}lung/loop/mustache
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
method=SCALE
chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)

for sample in ${samples_for_tissue[@]}
do
    mkdir -p ${data_path}${tissue}/loop/mustache/${sample}
    juicer_file=${data_path}${tissue}/juicer/${sample}.allValidPairs.hic
    for chr in ${chromosomes[@]}
    do
        mkdir -p ${data_path}${tissue}/loop/mustache/${sample}/${chr}
        python ${mustacheDir}mustache.py -f ${juicer_file} -ch ${chr} -r ${resolution} -pt 0.01 -norm ${method} -o ${data_path}${tissue}/loop/mustache/${sample}/${chr}/${sample}_${resolution}_${chr}_mustache_loop.tsv
    done
done