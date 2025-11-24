data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/WANG_cellular_aging_HiC/
HiC_Pro=/storage/zhangyanxiaoLab/suzhuojie/software/HiC-Pro_3.1.0/
samples=(G1 G2 DS1 DS2)
resolution=20000
mkdir -p ${data_path}ice_matrix
for sample in ${samples[@]}
do
    /storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/hicpro/bin/python ${HiC_Pro}scripts/ice --results_filename ${data_path}ice_matrix/${sample}_${resolution}_iced.matrix --filter_low_counts_perc 0.02 --filter_high_counts_perc 0 --max_iter 100 --eps 0.1 --remove-all-zeros-loci --output-bias 1 ${data_path}raw_matrix/${sample}_${resolution}.matrix &
done

HiC_Pro=/storage/zhangyanxiaoLab/suzhuojie/software/HiC-Pro_3.1.0/bin/utils/
chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
for sample in ${samples[@]}
do
    mkdir -p ${data_path}ice_matrix/${sample}
    for chr in ${chromosomes[@]} 
    do
        echo $sample $chr begin
        /storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/hicpro/bin/python ${HiC_Pro}split_sparse.py ${data_path}ice_matrix/${sample}_${resolution}_iced.matrix \
            -b ${data_path}raw_matrix/${sample}_${resolution}_abs.bed \
            -c $chr \
            -o ${data_path}ice_matrix/${sample}/${sample}_${resolution}_${chr}_iced.matrix
    done
done

for sample in ${samples[@]}
do 
    mkdir -p ${data_path}dense_matrix/${sample}
    for chr in ${chromosomes[@]}
    do
        echo $sample $chr begin
        /storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/hicpro/bin/python ${HiC_Pro}sparseToDense.py ${data_path}ice_matrix/${sample}/${sample}_${resolution}_${chr}_iced.matrix_${chr}.matrix \
            -b ${data_path}ice_matrix/${sample}/${sample}_${resolution}_${chr}_iced.matrix_${chr}_abs.bed \
            -i \
            -o ${data_path}dense_matrix/${sample}/${sample}_${resolution}_${chr}_dense.matrix
        gzip ${data_path}dense_matrix/${sample}/${sample}_${resolution}_${chr}_dense.matrix
    done
done