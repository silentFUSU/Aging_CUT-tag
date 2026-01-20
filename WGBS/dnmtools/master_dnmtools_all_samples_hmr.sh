tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
max_jobs=2
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do  
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/WGBS/dnmtools/dnmtools_all_samples_hmr.sh ${tissue} &
done