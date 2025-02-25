data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/WANG_cellular_aging_HiC/
scripts=/storage/zhangyanxiaoLab/suzhuojie/software/crane-nature-2015-master/scripts/
samples=(DS1 DS2 G1 G2)
resolution=20000
chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr22 chr20 chr21 chr22 chrX chrY)
mkdir -p ${data_path}TAD/
mkdir -p ${data_path}TAD/insulation_score

for sample in ${samples[@]}
do 
    mkdir -p ${data_path}TAD/insulation_score/${sample}
    cd ${data_path}TAD/insulation_score/${sample}
    for chr in ${chromosomes[@]}
    do
        echo $sample $chr begin
        perl ${scripts}/matrix2insulation.pl -i ${data_path}dense_matrix/${sample}/${sample}_${resolution}_${chr}_dense.matrix.gz \
            -is 500000 -ids 200000 -im mean -bmoe 3 -nt 0.1 -v
    done
done

for sample in ${samples[@]}
do 
    files=""
    for chr in {1..22} X Y
    do  
        files+=" ${data_path}TAD/insulation_score/${sample}/${sample}_${resolution}_chr${chr}_dense.is500001.ids200001.insulation.bedGraph"  
    done  
    awk 'FNR > 1' $files | awk -F'\t' '!(tolower($4) == "na")' > ${data_path}TAD/insulation_score/${sample}/${sample}_${resolution}_dense.is500001.ids200001.insulation.bedGraph
done

for sample in ${samples[@]}
do 
    files=""
    for chr in {1..22} X Y
    do  
        files+=" ${data_path}TAD/insulation_score/${sample}/${sample}_${resolution}_chr${chr}_dense.is500001.ids200001.insulation.boundaries.bed"  
    done  
    awk 'FNR > 1' $files  > ${data_path}TAD/insulation_score/${sample}/${sample}_${resolution}_dense.is500001.ids200001.insulation.boundaries.bed
done

for sample in ${samples[@]}
do 
    files=""
    for chr in {1..22} X Y
    do 
        files+=" ${data_path}TAD/insulation_score/${sample}/${sample}_${resolution}_chr${chr}_dense.is500001.ids200001.insulation"  
    done  
    awk 'FNR > 1' $files  | awk -F'\t' '!(tolower($9) == "na")' > ${data_path}TAD/insulation_score/${sample}/${sample}_${resolution}_dense.is500001.ids200001.insulation
done

for sample in ${samples[@]}
do 
    files=""
    for chr in {1..22} X Y
    do 
        files+=" ${data_path}TAD/insulation_score/${sample}/${sample}_${resolution}_chr${chr}_dense.is500001.ids200001.insulation.boundaries"  
    done  
    awk 'FNR > 1' $files  > ${data_path}TAD/insulation_score/${sample}/${sample}_${resolution}_dense.is500001.ids200001.insulation.boundaries
done