tissues=(lung liver brain CB kidney colon)
resolution=20000
code_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/HiC/
for tissue in ${tissues[@]}
do
    bash ${code_path}TAD/homer_findTADs.sh ${tissue} ${resolution} 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/${tissue}_homer_findTADs_loops_${resolution}.log
    bash ${code_path}Differential_analysis/homer_getDiffExpression.sh ${tissue} ${resolution} 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/${tissue}_homer_findTADs_loops_${resolution}.log
done
for tissue in ${tissues[@]}
do
    tail -n +2 ~/projects/Aging_CUT_Tag/data/samples/HiC/${tissue}/loop/homer/merged_${resolution}.loop.2D.bed | awk '{print $1, $2, $3, $4, $5, $6}' OFS='\t' > ~/projects/Aging_CUT_Tag/data/samples/HiC/${tissue}/loop/homer/merged_${resolution}.loop.2D.bedpe
done