tissues=(kidney ileum Hip)
for tissue in ${tissues[@]}
do  
    /usr/local/lib64/R/bin/Rscript ~/projects/Aging_CUT_Tag/code/WGBS/DMR/DSS_DMR.R ${tissue}
done
