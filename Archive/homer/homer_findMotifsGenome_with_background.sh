data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/ATAC/
# tissues=(bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
# tissues=(aorta BAT)
tissues=(iWAT)
conditions=(up down)
mkdir ${data_path}motif_bg
mkdir ${data_path}motif_bg/up
mkdir ${data_path}motif_bg/down
max_jobs=3
current_jobs() {  
    jobs -rp | wc -l  
}  
for tissue in ${tissues[@]}
do
    for condition in ${conditions[@]}
    do
        bed=${data_path}motif_bg/bed/${tissue}_ATAC_peaks_diff_${condition}_homer.bed
        if [ -s "$bed" ]; then  
            while [ $(current_jobs) -ge $max_jobs ]; do  
                sleep 1  
            done  
            echo ${tissue} ${condition} begin
            mkdir ${data_path}motif_bg/${condition}/${tissue}
            findMotifsGenome.pl ${bed} mm10 ${data_path}motif_bg/${condition}/${tissue} -size 200 -mask -bg ${data_path}motif_bg/bed/${tissue}_ATAC_peaks_diff_stable_homer.bed &
        fi
    done
done
wait
echo done