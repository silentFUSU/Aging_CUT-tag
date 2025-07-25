tissues=(brain CB kidney liver lung bonemarrow colon heart Hip mammarygland stomach thymus skin muscle)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
HiC_Pro=/storage/zhangyanxiaoLab/suzhuojie/software/HiC-Pro_3.1.0/
chromsize=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.chrom.sizes
result_path=/mnt/transposon2/zhangyanxiaoLab/suzhuojie/project/Aging_CUT_Tag/LMJ_Network_rawdata/HiC/
resolution=5000
# for tissue in ${tissues[@]}
# do
#     mkdir -p ${result_path}${tissue}
#     mkdir -p ${result_path}${tissue}/raw_matrix
#     search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
#     samples=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
#     for sample in ${samples[@]}
#     do
#         cat ${data_path}${tissue}/ValidPairs/${sample}.allValidPairs | ${HiC_Pro}scripts/build_matrix --matrix-format upper --binsize ${resolution} --chrsizes ${chromsize} --ifile /dev/stdin --oprefix ${result_path}${tissue}/raw_matrix/${sample}_${resolution}
#     done
# done

for tissue in ${tissues[@]}
do
    mkdir -p ${result_path}${tissue}
    mkdir -p ${result_path}${tissue}/ice_matrix 
    search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
    samples=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
    for sample in ${samples[@]}
    do 
        /storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/hicpro/bin/python ${HiC_Pro}scripts/ice --results_filename ${result_path}${tissue}/ice_matrix/${sample}_${resolution}_iced.matrix --filter_low_counts_perc 0.02 --filter_high_counts_perc 0 --max_iter 100 --eps 0.1 --remove-all-zeros-loci --output-bias 1 ${result_path}${tissue}/raw_matrix/${sample}_${resolution}.matrix
    done
done