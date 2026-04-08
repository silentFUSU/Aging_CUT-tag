# tissues=(brain lung kidney liver bonemarrow colon heart Hip mammarygland thymus)
tissues=muscle
for tissue in ${tissues[@]}
do
    snakemake -s /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/HiC/master_distance_normalize_snakefile \
        --cores 20 --config tissue="${tissue}" resolution="50000" 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/${tissue}_${resolution}_distance_normalize_snakefile.log 
done