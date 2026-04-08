tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
code_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/H3K27me3/
max_jobs=2
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    echo ${tissue} begin
    bash ${code_path}generate_fold_change_bw.sh $tissue 2>&1>>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/${tissue}_generate_fold_change_bw.log &
done
wait
echo all tissues done
