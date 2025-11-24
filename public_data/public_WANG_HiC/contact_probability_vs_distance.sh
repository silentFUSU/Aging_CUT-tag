data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/WANG_cellular_aging_HiC/
mkdir -p ${data_path}4DN_pairs
mkdir -p ${data_path}tmp
samples_for_tissue=(G1 G2 DS1 DS2)
resolution=20000
HiC_pro=/storage/zhangyanxiaoLab/suzhuojie/software/HiC-Pro_3.1.0/bin/utils/
chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
for sample in ${samples_for_tissue[@]}
do
    mkdir -p ${data_path}raw_matrix/${sample}
    for chr in ${chromosomes[@]} 
    do
    /storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/hicpro/bin/python ${HiC_pro}/split_sparse.py ${data_path}raw_matrix/${sample}_${resolution}.matrix \
        -b ${data_path}raw_matrix/${sample}_${resolution}_abs.bed \
        -c $chr \
        -o ${data_path}raw_matrix/${sample}/${sample}_${resolution}_${chr}_raw.matrix
    done
done

mkdir -p ${data_path}distance_contact/
res_k=$((${resolution}/1000))
for sample in ${samples_for_tissue[@]}
do
    files=""
    for chr in ${chromosomes[@]} 
    do
        files="$files ${data_path}raw_matrix/${sample}/${sample}_${resolution}_${chr}_raw.matrix_${chr}.matrix"
    done
    echo $files 
    awk -v res="$resolution" '{dist=($2-$1)*res;mat[dist]+=$3} END { for (dist in mat){print dist,mat[dist]} }' $files |sort -k1,1n >   ${data_path}distance_contact/$sample.dist.contacts.${res_k}k &
done
wait