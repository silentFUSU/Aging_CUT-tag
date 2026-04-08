data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/ATAC_peak_from_LMJ/
tissues=(skin aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas spleen stomach testis thymus tongue uterus mammarygland iWAT)
mkdir ${data_path}motif_without_bg
mkdir ${data_path}motif_without_bg/up
mkdir ${data_path}motif_without_bg/down
max_jobs=3
current_jobs() {  
    jobs -rp | wc -l  
}  

for tissue in ${tissues[@]}
do
    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    echo ${tissue} Up begin
    bed=${data_path}/up/${tissue}_Up_sorted_homer.bed
    mkdir ${data_path}motif_without_bg/up/${tissue}
    findMotifsGenome.pl ${bed} mm10 ${data_path}motif_without_bg/up/${tissue} -size 200 -mask 2>&1>~/projects/Aging_CUT_Tag/logs/${tissue}_ATAC_peak_MJ_up_findMotifsGenome_without_background.log &

    while [ $(current_jobs) -ge $max_jobs ]; do  
        sleep 1  
    done  
    echo ${tissue} Down begin
    bed=${data_path}/down/${tissue}_Down_sorted_homer.bed
    mkdir ${data_path}motif_without_bg/down/${tissue}
    findMotifsGenome.pl ${bed} mm10 ${data_path}motif_without_bg/down/${tissue} -size 200 -mask 2>&1>~/projects/Aging_CUT_Tag/logs/${tissue}_ATAC_peak_MJ_down_findMotifsGenome_without_background.log &
done
wait
echo done