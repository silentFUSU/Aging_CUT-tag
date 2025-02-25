data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=$1
res=$2
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table") 
mkdir -p ${data_path}${tissue}/TAD/
mkdir -p ${data_path}${tissue}/TAD/homer
mkdir -p ${data_path}${tissue}/loop/
mkdir -p ${data_path}${tissue}/loop/homer

for sample in ${samples_for_tissue[@]}
do
    mkdir -p ${data_path}${tissue}/TAD/homer/${sample}
    mkdir -p ${data_path}${tissue}/loop/homer/${sample}
    findTADsAndLoops.pl find ${data_path}${tissue}/compartment/homer_compartment/tagDir/${sample} -cpu 5 -res ${res} -genome mm10 \
        -badChr chr1_GL456210_random,chr1_GL456211_random,chr1_GL456212_random,chr1_GL456221_random,chr4_GL456216_random,chr4_GL456350_random,chr4_JH584292_random,chr4_JH584294_random,chr4_JH584295_random,chr5_GL456354_random,chr5_JH584296_random,chr5_JH584297_random,chr5_JH584298_random,chr5_JH584299_random,chr7_GL456219_random,chrM,chrUn_GL456239,chrUn_GL456359,chrUn_GL456360,chrUn_GL456366,chrUn_GL456367,chrUn_GL456368,chrUn_GL456370,chrUn_GL456372,chrUn_GL456378,chrUn_GL456379,chrUn_GL456381,chrUn_GL456382,chrUn_GL456383,chrUn_GL456385,chrUn_GL456387,chrUn_GL456389,chrUn_GL456390,chrUn_GL456392,chrUn_GL456393,chrUn_GL456394,chrUn_GL456396,chrUn_JH584304,chrX_GL456233_random,chrY_JH584300_random,chrY_JH584301_random,chrY_JH584303_random \
        -o ${data_path}${tissue}/TAD/homer/${sample}/${sample}_${res}
    mv ${data_path}${tissue}/TAD/homer/${sample}/${sample}_${res}.loop.2D.bed ${data_path}${tissue}/loop/homer/${sample}
done
