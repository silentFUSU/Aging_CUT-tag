tissues=(brain CB liver lung kidney colon)
code_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/HiC/
for tissue in ${tissues[@]}
do
    bash ${code_path}loop/HiCCUP_call_loop.sh ${tissue} 20000 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/${tissue}_20000_HiCCUP_call_loop.log &
    
done
for tissue in ${tissues[@]}
do
    bash ${code_path}loop/mustache_call_loop.sh ${tissue} 20kb 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/${tissue}_20kb_mustache_call_loop.log
done