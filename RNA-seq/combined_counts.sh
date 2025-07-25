data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
len=$(find ${data_path}*/counts/ -type f \( -name "LLX*.counts" -o -name "CKJ*.counts" -o -name "DYQ*.counts" \) ! -name "nodup.counts" | wc -l)
counts=$(ls ${data_path}*/counts/*.counts | grep -v "nodup.counts" | grep -E "^(${data_path}*/counts/(LLX|CKJ|DYQ)).*.counts")
paste ${counts} | cut -f 1-6,$(seq -s, 7 7 $((7*len))) | grep -v 'chrM' > ${data_path}combined-chrM.counts

tissues=(skin CB spleen heart bladder tongue uterus aorta thymus stomach Hip FC BAT iWAT muscle bonemarrow lung kidney liver testis colon cecum ileum jejunum ovary mammarygland)
for tissue in ${tissues[@]}
do
    data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
    len=$(find ${data_path}${tissue}/counts/ -type f -name *.counts | grep -v nodup.counts | wc -l)  
    counts=$(ls ${data_path}${tissue}/counts/*.counts | grep -v nodup.counts )
    paste ${counts} | cut -f 1-6,$(seq -s, 7 7 $((7*len))) | grep -v 'chrM' > ${data_path}${tissue}/combined-chrM.counts
done
