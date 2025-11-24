data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/WANG_cellular_aging_HiC/
HiC_Pro=/storage/zhangyanxiaoLab/suzhuojie/software/HiC-Pro_3.1.0/bin/utils/
samples=(G1 G2 DS1 DS2)
resolution=200000
mkdir -p ${data_path}dense_matrix
for sample in ${samples[@]}
do
    mkdir -p ${data_path}dense_matrix/${sample}
    mkdir -p ${data_path}dense_matrix/${sample}/raw
    echo $sample begin
    cd ${data_path}dense_matrix/${sample}/raw/
    /storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/hicpro/bin/python ${HiC_Pro}sparseToDense.py ${data_path}raw_matrix/${sample}_${resolution}.matrix \
        -b ${data_path}raw_matrix/${sample}_${resolution}_abs.bed \
        -c
done