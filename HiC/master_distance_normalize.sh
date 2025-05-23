tissues=(brain CB kidney liver lung bonemarrow colon heart Hip mammarygland stomach thymus)
max_jobs=2
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/HiC/distance_normalize.sh ${tissue} 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/${tissue}_HiC_ob_ex_make.log &
done