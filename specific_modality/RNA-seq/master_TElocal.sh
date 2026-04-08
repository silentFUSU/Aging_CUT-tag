data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
species=mm10
max_jobs=5
current_jobs() {  
    jobs -rp | wc -l  
}  

for tissue in ${tissues[@]}
do
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    echo ${tissue} TElocal begin
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/RNA-seq/TElocal.sh  -i ${data_path} -g ${species} -t ${tissue} &
done
wait