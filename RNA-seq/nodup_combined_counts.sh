data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
len=$(find ${data_path}*/counts/ -type f -name "*.nodup.counts" | wc -l)  
counts=$(ls ${data_path}*/counts/*.nodup.counts )
paste ${counts} | cut -f 1-6,$(seq -s, 7 7 $((7*len))) | grep -v 'chrM' > ${data_path}combined-chrM.nodup.counts