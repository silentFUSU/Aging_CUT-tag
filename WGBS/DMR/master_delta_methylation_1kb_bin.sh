# tissues=(liver lung mammarygland kidney ileum Hip skin bonemarrow jejunum colon ovary CB BAT thymus testis heart muscle stomach bladder aorta tongue)
tissues=(bladder aorta tongue)
max_jobs=2
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do  
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    bash ~/projects/Aging_CUT_Tag/code/WGBS/DMR/delta_methylation_1kb_bin.sh ${tissue} &
done
wait

# for tissue in ${tissues[@]}
# do  
#     while [ $(current_jobs) -ge $max_jobs ]; do  
#         sleep 1  
#     done  
#     bash ~/projects/Aging_CUT_Tag/code/WGBS/DMR/WGBS_1kb_delta_in_DEG_body_plotprofile.sh ${tissue} &
# done
# wait