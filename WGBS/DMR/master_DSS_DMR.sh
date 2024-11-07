# tissues=()
# tissues=(skin mammarygland BAT thymus testis liver lung  kidney ileum Hip bonemarrow jejunum colon ovary CB muscle heart stomach baldder aorta tongue)
tissues=(bladder)
max_jobs=2
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