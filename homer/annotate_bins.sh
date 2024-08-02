data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/ChromHMM/15_until_skin/state_transfer/E9_to_E1/
ref_path=/storage/zhangyanxiaoLab/share/
mkdir ${data_path}annotate/
tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum ileum bonemarrow ileum heart thymus stomach skin bladder tongue aorta)
for tissue in ${tissues[@]}
do
    annotatePeaks.pl ${data_path}bed/${tissue}_homer.bed mm10 -gtf ${ref_path}gtf/mm10* > ${data_path}annotate/${tissue}_homer.txt
done