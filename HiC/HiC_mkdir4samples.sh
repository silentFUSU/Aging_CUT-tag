tissue=$1
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
mkdir ${data_path}${tissue}
files=(juicer ValidPairs)
for file in ${files[@]}
do
    mkdir ${data_path}${tissue}/${file}
done