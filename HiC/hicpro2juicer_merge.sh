data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240925_WJH_Lung_HiC/result/hic_results/data/
HiC_Pro=/storage/zhangyanxiaoLab/suzhuojie/software/HiC-Pro_3.1.0/
chromsize=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.chrom.sizes
restriction=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/restriction_sites/mm10_dpnii.bed
juicer_tools=/storage/zhangyanxiaoLab/suzhuojie/software/juicer/scripts/common/juicer_tools.jar
samples=$(find ${data_path} -mindepth 1 -maxdepth 1 -type d -printf '%f\n')  
for sample in ${samples[@]}
do
    bash ${HiC_Pro}bin/utils/hicpro2juicebox.sh -i ${data_path}${sample}/${sample}.allValidPairs -g ${chromsize} -j ${juicer_tools} -r ${restriction} -o ${data_path} -t ${data_path}tmp &
done
wait