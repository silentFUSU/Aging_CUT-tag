data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=lung
resolution=1mb
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
ref=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.fa
mkdir -p ${data_path}${tissue}/compartment
for sample in ${samples_for_tissue[@]}
do
    juicer_hic=${data_path}${tissue}/juicer/${sample}.allValidPairs.hic
    fanc compartments -g ${ref} \
                  -d ${data_path}${tissue}/compartment/fanc_${sample}_${resolution}.compartment.bed \
                  -i 2 -v ${data_path}${tissue}/compartment/fanc_${sample}_${resolution}.ev.txt \
                  ${juicer_hic}@${resolution}  ${data_path}${tissue}/compartment/fanc_${sample}_${resolution}.ab &
done
wait 
echo done