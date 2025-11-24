tissues=(BAT mammarygland CB lung kidney aorta brain spleen thymus skin bladder bonemarrow Hip heart muscle jejunum uterus ovary liver tongue cecum colon testis stomach pancreas iWAT ileum)
max_jobs=1
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    echo ${tissue} motif begin
    python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/WGBS/DMR/snapatac2_DMR_motif.py --tissue ${tissue} &
done