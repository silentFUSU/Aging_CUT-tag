# tissues=(skin mammarygland BAT thymus testis liver lung  kidney ileum Hip bonemarrow jejunum colon ovary CB muscle heart stomach bladder aorta tongue)
tissues=(BAT mammarygland testis thymus)
code_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/WGBS/DMR/
max_jobs=2
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    # bash ${code_path}WGBS_change_in_histone_change_plotprofile.sh ${tissue}
    # bash ${code_path}WGBS_change_in_DEG_TSS_plotprofile.sh ${tissue}
    bash ${code_path}WGBS_change_in_DEG_body_plotprofile.sh ${tissue} &
done