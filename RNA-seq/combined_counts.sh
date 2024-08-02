data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
len=$(ls ${data_path}*/counts/*.counts | wc -l)  
counts=$(ls ${data_path}*/counts/*.counts)
paste ${counts} | cut -f 1-6,$(seq -s, 7 7 $((7*len))) | grep -v 'chrM' > ${data_path}combined-chrM.counts