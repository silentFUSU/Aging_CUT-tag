tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
max_jobs=2
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do  
    echo $tissue begin
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    bash ~/projects/Aging_CUT_Tag/code/WGBS/DMR/delta_methylation_1kb_bin.sh ${tissue} &
done
wait
echo all done
# for tissue in ${tissues[@]}
# do  
#     while [ $(current_jobs) -ge $max_jobs ]; do  
#         sleep 1  
#     done  
#     bash ~/projects/Aging_CUT_Tag/code/WGBS/DMR/WGBS_1kb_delta_in_DEG_body_plotprofile.sh ${tissue} &
# done
# wait