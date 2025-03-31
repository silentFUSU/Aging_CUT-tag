data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=$1
chromsize=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.normal.chrom.sizes
pairix=/storage/zhangyanxiaoLab/suzhuojie/software/pairix/bin/
mkdir -p ${data_path}${tissue}/4DN_pairs
mkdir -p ${data_path}${tissue}/tmp
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
resolution=$2
HiC_pro=/storage/zhangyanxiaoLab/suzhuojie/software/HiC-Pro_3.1.0/bin/utils/
chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
for sample in ${samples_for_tissue[@]}
do
    mkdir -p ${data_path}${tissue}/raw_matrix/${sample}
    for chr in ${chromosomes[@]} 
    do
    /storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/hicpro/bin/python ${HiC_pro}/split_sparse.py ${data_path}${tissue}/raw_matrix/${sample}_${resolution}.matrix \
        -b ${data_path}${tissue}/raw_matrix/${sample}_${resolution}_abs.bed \
        -c $chr \
        -o ${data_path}${tissue}/raw_matrix/${sample}/${sample}_${resolution}_${chr}_raw.matrix
    done
done

mkdir -p ${data_path}${tissue}/distance_contact/
for sample in ${samples_for_tissue[@]}
do
    files=""
    for chr in ${chromosomes[@]} 
    do
        files="$files ${data_path}${tissue}/raw_matrix/${sample}/${sample}_${resolution}_${chr}_raw.matrix_${chr}.matrix"
    done
    echo $files 
    awk '{dist=($2-$1)*10000;mat[dist]+=$3} END { for (dist in mat){print dist,mat[dist]} }' $files |sort -k1,1n >   ${data_path}${tissue}/distance_contact/$sample.dist.contacts.10k &
done
wait