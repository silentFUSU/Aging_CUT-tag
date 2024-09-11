tissue=lung
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/WGBS/
mkdir ${data_path}${tissue}
files=(bdg bw)
for file in ${files[@]}
do
    mkdir ${data_path}${tissue}/${file}
done