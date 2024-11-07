data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=lung
resolution=10000
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
samples_for_tissue=(WJH-103-Lung WJH-106-Lung WJH-109-Lung)
mkdir -p ${data_path}${tissue}/cool
source /storage/zhangyanxiaoLab/suzhuojie/miniconda3/etc/profile.d/conda.sh
conda activate hicexplorer

for sample in ${samples_for_tissue[@]}
do
    if [ ! -f "${data_path}${tissue}/cool/${sample}.allValidPairs.${resolution}.cool" ]; then 
        hic2cool convert ${data_path}${tissue}/juicer/${sample}.allValidPairs.hic ${data_path}${tissue}/cool/${sample}.allValidPairs.${resolution}.cool -r ${resolution}  
        cooler balance ${data_path}${tissue}/cool/${sample}.allValidPairs.${resolution}.cool
    fi
done

source /storage/zhangyanxiaoLab/suzhuojie/miniconda3/etc/profile.d/conda.sh
conda activate domaincaller
mkdir -p ${data_path}${tissue}/TAD/domaincaller
for sample in ${samples_for_tissue[@]}
do
    if [ -f "${data_path}${tissue}/cool/${sample}.allValidPairs.${resolution}.cool" ]; then 
        domaincaller --uri ${data_path}${tissue}/cool/${sample}.allValidPairs.${resolution}.cool -O ${data_path}${tissue}/TAD/domaincaller/${sample}.allValidPairs.${resolution}.output \
            -D ${data_path}${tissue}/TAD/domaincaller/${sample}.allValidPairs.${resolution}.DI -p 1 \
            --exclude chrM chrL chrM chr1_GL456210_random chr1_GL456211_random chr1_GL456212_random chr1_GL456213_random chr1_GL456221_random chr4_GL456216_random chr4_GL456350_random chr4_JH584292_random chr4_JH584293_random chr4_JH584294_random chr4_JH584295_random chr5_GL456354_random chr5_JH584296_random chr5_JH584297_random chr5_JH584298_random chr5_JH584299_random chr7_GL456219_random chrX_GL456233_random chrY_JH584300_random chrY_JH584301_random chrY_JH584302_random chrY_JH584303_random chrUn_GL456239 chrUn_GL456359 chrUn_GL456360 chrUn_GL456366 chrUn_GL456367 chrUn_GL456368 chrUn_GL456370 chrUn_GL456372 chrUn_GL456378 chrUn_GL456379 chrUn_GL456381 chrUn_GL456382 chrUn_GL456383 chrUn_GL456385 chrUn_GL456387 chrUn_GL456389 chrUn_GL456390 chrUn_GL456392 chrUn_GL456393 chrUn_GL456394 chrUn_GL456396 chrUn_JH584304 \
            --logFile /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/${sample}_domaincaller.log 
        awk -F',' '{print $1 "\t" $2 "\t" $3}' ${data_path}${tissue}/TAD/domaincaller/${sample}.allValidPairs.${resolution}.output > ${data_path}${tissue}/TAD/domaincaller/${sample}.allValidPairs.${resolution}.bed
    fi
done

echo all done