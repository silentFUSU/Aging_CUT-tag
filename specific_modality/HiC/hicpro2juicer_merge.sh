tissue=$1
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
HiC_Pro=/storage/zhangyanxiaoLab/suzhuojie/software/HiC-Pro_3.1.0/
chromsize=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.chrom.sizes
samples=$(ls ${data_path}${tissue}/ValidPairs/*.allValidPairs |  sed 's|.*/||; s|\.allValidPairs$||') 
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
juicer_tools=/storage/zhangyanxiaoLab/suzhuojie/software/juicer/scripts/common/juicer_tools.jar
mkdir -p ${data_path}${tissue}/tmp
for sample in ${samples[@]}
do
    if [ ! -f "${data_path}${tissue}/juicer/${sample}.allValidPairs.hic" ]; then
        echo ${sample}
        bash ${HiC_Pro}bin/utils/hicpro2juicebox.sh -i ${data_path}${tissue}/ValidPairs/${sample}.allValidPairs -g ${chromsize} -j ${juicer_tools} -o ${data_path}${tissue}/juicer/ -t ${data_path}${tissue}/tmp &
    fi  
done
wait
echo done

# data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20241205_WJH_Heart_HiC/result/hic_results/data/
# HiC_Pro=/storage/zhangyanxiaoLab/suzhuojie/software/HiC-Pro_3.1.0/
# chromsize=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.chrom.sizes
# samples=$(ls ${data_path}/*/*.allValidPairs |  sed 's|.*/||; s|\.allValidPairs$||') 
# juicer_tools=/storage/zhangyanxiaoLab/suzhuojie/software/juicer/scripts/common/juicer_tools.jar
# mkdir -p ${data_path}/tmp
# for sample in ${samples[@]}
# do
#     bash ${HiC_Pro}bin/utils/hicpro2juicebox.sh -i ${data_path}${sample}/${sample}.allValidPairs -g ${chromsize} -j ${juicer_tools} -o ${data_path} -t ${data_path}/tmp &  
# done