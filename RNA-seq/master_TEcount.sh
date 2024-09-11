data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
# tissues=(jejunum ileum cecum colon testis liver kidney lung bonemarrow muscle iWAT)
# tissues=(pancreas skin CB spleen heart bladder tongue uterus aorta thymus stomach Hip FC BAT iWAT muscle bonemarrow lung kidney liver testis colon cecum ileum jejunum)
tissues=(mammarygland)
species=mm10
max_jobs=7
current_jobs() {  
    jobs -rp | wc -l  
}  

for tissue in ${tissues[@]}
do
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    echo ${tissue} TEcount begin
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/RNA-seq/TEcount.sh  -i ${data_path} -g ${species} -t ${tissue} &
done