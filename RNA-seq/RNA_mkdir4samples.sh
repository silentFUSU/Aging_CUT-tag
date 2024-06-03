tissue=$1
dirs=(bam counts bw)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/

mkdir ${data_path}${tissue}

for dir in ${dirs[@]}
do
    mkdir ${data_path}${tissue}/${dir}
done 

