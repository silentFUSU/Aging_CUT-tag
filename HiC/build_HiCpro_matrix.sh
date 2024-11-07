tissue=lung
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
HiC_Pro=/storage/zhangyanxiaoLab/suzhuojie/software/HiC-Pro_3.1.0/
chromsize=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.chrom.sizes
samples=$(ls ${data_path}${tissue}/ValidPairs/WJH*.allValidPairs | sed 's|.*/||; s|\.allValidPairs$||' ) 
# resolution=10000
resolution=50000
mkdir -p ${data_path}${tissue}/raw_matrix
# samples=(WJH-103-Lung WJH-106-Lung WJH-109-Lung)
for sample in ${samples[@]}
do
    cat ${data_path}${tissue}/ValidPairs/${sample}.allValidPairs | ${HiC_Pro}scripts/build_matrix --matrix-format upper --binsize ${resolution} --chrsizes ${chromsize} --ifile /dev/stdin --oprefix ${data_path}${tissue}/raw_matrix/${sample}_${resolution}
done

mkdir -p ${data_path}${tissue}/ice_matrix
for sample in ${samples[@]}
do
    /storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/hicpro/bin/python ${HiC_Pro}scripts/ice --results_filename ${data_path}${tissue}/ice_matrix/${sample}_${resolution}_iced.matrix --filter_low_counts_perc 0.02 --filter_high_counts_perc 0 --max_iter 100 --eps 0.1 --remove-all-zeros-loci --output-bias 1 ${data_path}${tissue}/raw_matrix/${sample}_${resolution}.matrix
done