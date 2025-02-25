tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
max_jobs=5
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do  

    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    echo $tissue macs2 begin
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/call_peak/macs2_age_split_bam_callpeaks.sh $tissue 2>&1>~/projects/Aging_CUT_Tag/logs/${tissue}_macs2_age_split_bam_callpeaks.log &

    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    echo $tissue sicer2 begin
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/call_peak/sicer2_age_split_bam_callpeaks.sh $tissue 2>&1>~/projects/Aging_CUT_Tag/logs/${tissue}_sicer2_age_split_bam_callpeaks.log &
done
wait 
echo ${tissues[@]} all done

