source /storage/zhangyanxiaoLab/suzhuojie/miniconda3/etc/profile.d/conda.sh
conda activate hicexplorer
tissue=brain
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
samples=$(ls ${data_path}${tissue}/ValidPairs/{WJH,DYQ}*.allValidPairs | sed 's|.*/||; s|\.allValidPairs$||' ) 
resolution=50000
mkdir -p ${data_path}${tissue}/cool/
for sample in ${samples[@]}
do
    hic2cool convert ${data_path}${tissue}/juicer/${sample}.allValidPairs.hic ${data_path}${tissue}/cool/${sample}_${resolution}.cool -r ${resolution} &
done
wait 
echo all done