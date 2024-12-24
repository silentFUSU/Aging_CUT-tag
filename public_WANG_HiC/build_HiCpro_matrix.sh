data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/WANG_cellular_aging_HiC/
HiC_Pro=/storage/zhangyanxiaoLab/suzhuojie/software/HiC-Pro_3.1.0/
chromsize=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/hg38/hg38.chrom.sizes
resolution=200000
samples=(G1 G2 DS1 DS2)
mkdir -p ${data_path}/raw_matrix
for sample in ${samples[@]}
do
    cat ${data_path}ValidPairs/${sample}.allValidPairs | ${HiC_Pro}scripts/build_matrix --matrix-format upper --binsize ${resolution} --chrsizes ${chromsize} --ifile /dev/stdin --oprefix ${data_path}raw_matrix/${sample}_${resolution} &
done