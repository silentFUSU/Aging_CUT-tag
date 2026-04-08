tissue=$1
data_path=/mnt/transposon2/zhangyanxiaoLab/suzhuojie/project/Aging_CUT_Tag/samples/HiC/
mkdir ${data_path}${tissue}
files=(juicer ValidPairs)
for file in ${files[@]}
do
    mkdir ${data_path}${tissue}/${file}
done
