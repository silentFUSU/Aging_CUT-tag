tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
max_jobs=4
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do  
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    /usr/local/lib64/R/bin/Rscript ~/projects/Aging_CUT_Tag/code/WGBS/DMR/DSS_DMR.R ${tissue} &
done
wait
echo done