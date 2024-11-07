data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=lung
scripts=/storage/zhangyanxiaoLab/suzhuojie/software/crane-nature-2015-master/scripts/
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table") 
resolution=10000
chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
mkdir -p ${data_path}${tissue}/TAD/insulation_score

for sample in ${samples_for_tissue[@]}
do 
    mkdir -p ${data_path}${tissue}/TAD/insulation_score/${sample}
    cd ${data_path}${tissue}/TAD/insulation_score/${sample}
    for chr in ${chromosomes[@]}
    do
        echo $sample $chr begin
        perl ${scripts}/matrix2insulation.pl -i ${data_path}${tissue}/dense_matrix/${sample}/${sample}_${resolution}_${chr}_dense.matrix.gz \
            -is 500000 -ids 200000 -im mean -bmoe 3 -nt 0.1 -v
    done
done