tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
antibodys=(H3K27ac H3K4me3 H3K4me1 H3K36me3)
code_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/diff_peak/
max_jobs=4
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        while [ $(current_jobs) -ge $max_jobs ]; do  
            sleep 1  
        done  
        echo ${tissue} ${antibody} begin
        bash ${code_path}generate_fold_change_bw.sh $tissue $antibody 2>&1>>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/${tissue}_${antibody}_generate_fold_change_bw.log &
    done
done
wait
echo all tissues done
