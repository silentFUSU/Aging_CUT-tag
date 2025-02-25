source /storage/zhangyanxiaoLab/suzhuojie/miniconda3/etc/profile.d/conda.sh
conda activate cooltools
tissue=$1
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
samples=$(ls ${data_path}${tissue}/ValidPairs/{WJH,DYQ}*.allValidPairs | sed 's|.*/||; s|\.allValidPairs$||' ) 
resolution=50000
for sample in ${samples[@]}
do
    cooltools rearrange --view /storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.chrom.withoutchrY.bed \
        --new-chrom-col new_chrom \
        ${data_path}${tissue}/cool/${sample}_${resolution}.cool \
        ${data_path}${tissue}/cool/${sample}_${resolution}_without_chrY.cool
done
wait
echo all done