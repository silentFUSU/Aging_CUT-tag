data_path=/storage/zhangyanxiaoLab/fastq/2024/2024-02-02-HaploX-SZJ/
out_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/20240205_SZJ/
mkdir ${out_path}
mkdir ${out_path}raw_data
file_list=$(ls ${data_path}/SZJ{161,162,175,176,177,183,189,190,191}* | grep -o 'SZJ[0-9]\+' | sort -u) 
for file in ${file_list[@]}
do
mkdir ${out_path}raw_data/${file}
ln -s ${data_path}${file}* ${out_path}raw_data/${file} 
done