tissues=(brain CB kidney liver lung bonemarrow colon heart Hip mammarygland stomach thymus)
code_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/HiC/
max_jobs=3
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    echo ${tissue} HiCCUP begin
    bash ${code_path}loop/HiCCUP_call_loop.sh ${tissue} 25000 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/${tissue}_25000_HiCCUP_call_loop.log &
done
for tissue in ${tissues[@]}
do
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    echo ${tissue} mustache begin
    bash ${code_path}loop/mustache_call_loop.sh ${tissue} 25kb 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/${tissue}_25kb_mustache_call_loop.log &
done
wait
echo all done