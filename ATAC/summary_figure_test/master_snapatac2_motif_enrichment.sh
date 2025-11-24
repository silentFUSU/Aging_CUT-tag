tissues=(BAT mammarygland CB lung kidney aorta brain spleen thymus skin bladder bonemarrow Hip heart muscle jejunum uterus ovary liver tongue cecum colon testis stomach pancreas iWAT ileum)
conditions=(up down)
max_jobs=1
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do
    for condition in ${conditions[@]}
    do
        while [ $(current_jobs) -ge $max_jobs ]; do  
            sleep 1  
        done  
        echo ${tissue} ${condition} motif begin
        python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/ATAC/summary_figure_test/snapatac2_motif_enrichment.py --condition ${condition} --tissue ${tissue} &
    done
done