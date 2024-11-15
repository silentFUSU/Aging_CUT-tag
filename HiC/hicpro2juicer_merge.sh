tissue=$1
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
HiC_Pro=/storage/zhangyanxiaoLab/suzhuojie/software/HiC-Pro_3.1.0/
chromsize=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.chrom.sizes
samples=$(ls ${data_path}${tissue}/ValidPairs/*.allValidPairs |  sed 's|.*/||; s|\.allValidPairs$||') 
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
juicer_tools=/storage/zhangyanxiaoLab/suzhuojie/software/juicer/scripts/common/
mkdir -p ${data_path}${tissue}/tmp
for sample in ${samples[@]}
do
    if [ ! -f "${data_path}${tissue}/juicer/${sample}.allValidPairs.hic" ]; then
        echo ${sample}
        bash ${HiC_Pro}bin/utils/hicpro2juicebox.sh -i ${data_path}${tissue}/ValidPairs/${sample}.allValidPairs -g ${chromsize} -j ${juicer_tools} -o ${data_path}${tissue}/juicer/ -t ${data_path}${tissue}/tmp &
    fi  
done
wait
echo done


# juiceDir="/storage/zhangyanxiaoLab/suzhuojie/software/juicer/"
# for sample in ${samples[@]}
# do
#     cp ${data_path}${tissue}/juicer/${sample}.allValidPairs.hic ${data_path}${tissue}/juicer/${sample}.allValidPairs.KR.hic 
#     ${juiceDir}/scripts/common/juicer_tools addNorm ${data_path}${tissue}/juicer/${sample}.allValidPairs.KR.hic -w 1000 -j 5 -k KR
# done


# cleaned_file=$(mktemp) 
# cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
# young=$(awk -v tissue="$tissue" -F, '$1 == tissue && $5 == "3M" {print $3}' "$cleaned_file")  
# old=$(awk -v tissue="$tissue" -F, '$1 == tissue && $5 == "24M" {print $3}' "$cleaned_file")  

# echo "$young" | xargs -n 1 -I {} cat "${data_path}${tissue}/ValidPairs/{}.allValidPairs" > "${data_path}${tissue}/ValidPairs/young_combined.allValidPairs" &
# echo "$old" | xargs -n 1 -I {} cat "${data_path}${tissue}/ValidPairs/{}.allValidPairs" > "${data_path}${tissue}/ValidPairs/old_combined.allValidPairs" &
# wait
# bash ${HiC_Pro}bin/utils/hicpro2juicebox.sh -i ${data_path}${tissue}/ValidPairs/young_combined.allValidPairs -g ${chromsize} -j ${juicer_tools} -o ${data_path}${tissue}/juicer/ -t ${data_path}${tissue}/tmp &
# bash ${HiC_Pro}bin/utils/hicpro2juicebox.sh -i ${data_path}${tissue}/ValidPairs/old_combined.allValidPairs -g ${chromsize} -j ${juicer_tools} -o ${data_path}${tissue}/juicer/ -t ${data_path}${tissue}/tmp &

# wait
# echo all done 