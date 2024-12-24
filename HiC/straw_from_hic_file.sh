data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/yanxiao_2019NG_HiC/
result_path=${data_path}juicer_matrix
mkdir -p ${result_path}
samples=(D00_Rep1 D00_Rep2 D80_Rep1 D80_Rep2)

for sample in ${samples[@]}
do
    straw NONE ${data_path}juicer/${sample}.hic 
done