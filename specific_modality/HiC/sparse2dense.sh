tissue=$1
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
HiC_Pro=/storage/zhangyanxiaoLab/suzhuojie/software/HiC-Pro_3.1.0/bin/utils/
chromsize=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.chrom.sizes
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table") 
resolution=$2
mkdir -p ${data_path}${tissue}/dense_matrix
chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
for sample in ${samples[@]}
do
    mkdir -p ${data_path}${tissue}/ice_matrix/${sample}
    for chr in ${chromosomes[@]} 
    do
        echo $sample $chr begin
        /storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/hicpro/bin/python ${HiC_Pro}split_sparse.py ${data_path}${tissue}/ice_matrix/${sample}_${resolution}_iced.matrix \
            -b ${data_path}${tissue}/raw_matrix/${sample}_${resolution}_abs.bed \
            -c $chr \
            -o ${data_path}${tissue}/ice_matrix/${sample}/${sample}_${resolution}_${chr}_iced.matrix
    done
done

for sample in ${samples[@]}
do 
    mkdir -p ${data_path}${tissue}/dense_matrix/${sample}
    for chr in ${chromosomes[@]}
    do
        echo $sample $chr begin
        /storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/hicpro/bin/python ${HiC_Pro}sparseToDense.py ${data_path}${tissue}/ice_matrix/${sample}/${sample}_${resolution}_${chr}_iced.matrix_${chr}.matrix \
            -b ${data_path}${tissue}/ice_matrix/${sample}/${sample}_${resolution}_${chr}_iced.matrix_${chr}_abs.bed \
            -i \
            -o ${data_path}${tissue}/dense_matrix/${sample}/${sample}_${resolution}_${chr}_dense.matrix
        gzip ${data_path}${tissue}/dense_matrix/${sample}/${sample}_${resolution}_${chr}_dense.matrix
    done
done
echo all done

# for sample in ${samples[@]}
# do 
#     mkdir -p ${data_path}${tissue}/dense_matrix/${sample}/raw
#     echo $sample begin
#     cd ${data_path}${tissue}/dense_matrix/${sample}/raw/
#     /storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/hicpro/bin/python ${HiC_Pro}sparseToDense.py ${data_path}${tissue}/raw_matrix/${sample}_${resolution}.matrix \
#         -b ${data_path}${tissue}/raw_matrix/${sample}_${resolution}_abs.bed \
#         -c 
#     # gzip ${data_path}${tissue}/dense_matrix/${sample}/${sample}_${resolution}_raw_dense.matrix
# done
# echo all done
