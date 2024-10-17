# tissues=(brain liver testis colon kidney lung spleen muscle Hip cecum bonemarrow heart thymus stomach skin aorta tongue bladder CB jejunum uterus ovary ileum pancreas)
# tissues=(BAT)
tissues=(iWAT)
for tissue in ${tissues[@]}
do  
    echo $tissue
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/call_peak/macs2_age_split_bam_callpeaks.sh $tissue &
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/call_peak/sicer2_age_split_bam_callpeaks.sh $tissue &
done
wait 
echo ${tissues[@]} all done

