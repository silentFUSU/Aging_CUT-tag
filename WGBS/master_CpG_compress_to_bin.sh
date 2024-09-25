tissues=(liver lung mammarygland kidney ileum Hip)
max_jobs=2
current_jobs() {  
    jobs -rp | wc -l  
}  

for tissue in ${tissues[@]}
do
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    /usr/local/lib64/R/bin/Rscript ~/projects/Aging_CUT_Tag/code/WGBS/CpG_compress_to_bin.R ${tissue} 1kb  &
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    /usr/local/lib64/R/bin/Rscript ~/projects/Aging_CUT_Tag/code/WGBS/CpG_compress_to_bin.R ${tissue} 10kb &
done
wait