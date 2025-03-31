raw_data=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data_transposon2/20250303_WJH_HiC/
data_path=/mnt/transposon2/zhangyanxiaoLab/suzhuojie/project/Aging_CUT_Tag/samples/HiC/
search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv
samples=$(find ${raw_data}result/hic_results/data/ -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | grep -v tmp)  
cleaned_file=$(mktemp)  
cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
tissue_array=() 
while IFS= read -r prefix; do  
    tissue=$(awk -v prefix="$prefix" -F, '$3 == prefix {print $1}' "$cleaned_file")  
    if [ -n "$tissue" ]; then  
        tissue_array+=("$tissue")  
    fi  
done <<< "$samples" 
unique_tissue_array=($(printf "%s\n" "${tissue_array[@]}" | sort -u))  

for tissue in ${unique_tissue_array[@]}
do
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/HiC/HiC_mkdir4samples.sh ${tissue}
    samples_for_tissue=$(awk -F',' -v t="$tissue" 'NR > 1 && ($1 == t) {print $3}' "$search_table")  
    for sample_for_tissue in ${samples_for_tissue[@]}
    do
        if ls ${raw_data}result/hic_results/data/${sample_for_tissue}/${sample_for_tissue}*.allValidPairs 1> /dev/null 2>&1; then  
            ln -s ${raw_data}result/hic_results/data/${sample_for_tissue}/${sample_for_tissue}*.allValidPairs ${data_path}${tissue}/ValidPairs/
        fi
        if ls ${raw_data}result/hic_results/data/${sample_for_tissue}*.allValidPairs.hic 1> /dev/null 2>&1; then  
            ln -s ${raw_data}result/hic_results/data/${sample_for_tissue}*.allValidPairs.hic ${data_path}${tissue}/juicer/
        fi
    done
done

for tissue in ${unique_tissue_array[@]}
do  
    ln -s ${data_path}${tissue} /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
done

code_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/HiC/
log_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/
#### valid pairs to juicer
for tissue in ${unique_tissue_array[@]}
do
    bash ${code_path}hicpro2juicer_merge.sh $tissue 2>&1>${log_path}${tissue}_hicpro2juicer_merge.log &
done

### juicer to cooler
resolution=50000
for tissue in ${unique_tissue_array[@]}
do
    nohup bash ${code_path}juicer2cool.sh ${tissue} ${resolution} 2>&1>${log_path}${tissue}_juicer2cool_${resolution}.log &
done

### homer compartment
for tissue in ${unique_tissue_array[@]}
do
    bash ${code_path}compartment/homer_compartment.sh $tissue 2>&1>${log_path}${tissue}_homer_compartment.log 
done

### saddle plot 
resolution=50000
for tissue in ${unique_tissue_array[@]}
do
    nohup python ${code_path}compartment/compartment_saddleplot.py -t ${tissue} -r ${resolution} 2>&1>${log_path}${tissue}_compartment_saddleplot_homer_${resolution}.log &
done

###make matrix 
resolutions=(10000 20000 200000)
for tissue in ${unique_tissue_array[@]}
do
    for resolution in ${resolutions[@]}
    do
        bash ${code_path}build_HiCpro_matrix.sh ${tissue} ${resolution} 2>&1>${log_path}${tissue}_build_HiCpro_matrix_${resolution}.log
    done
done

### sparse2dense
resolutions=(10000 20000)
for tissue in ${unique_tissue_array[@]}
do  
    for resolution in ${resolutions[@]}
    do
       nohup bash ${code_path}sparse2dense.sh ${tissue} ${resolution} 2>&1>${log_path}${tissue}_sparse2dense_${resolution}.log
    done
done

resolutions=(200000)
for tissue in ${unique_tissue_array[@]}
do  
    for resolution in ${resolutions[@]}
    do
       nohup bash ${code_path}sparse2dense_binomial.sh ${tissue} ${resolution} 2>&1>${log_path}${tissue}_sparse2dense_binomial_${resolution}.log
    done
done

### contact probability vs distance
resolution=10000
for tissue in ${unique_tissue_array[@]}
do
    nohup bash ${code_path}quality_control/contact_probability_vs_distance.sh ${tissue} ${resolution} 2>&1>${log_path}${tissue}_contact_probability_vs_distance_${resolution}.log &
done

### binomial differential analysis
resolutions=(200000)
for tissue in ${unique_tissue_array[@]}
do  
    for resolution in ${resolutions[@]}
    do
       nohup bash ${code_path}Differential_analysis/Wang_diff_analysis_master.sh ${tissue} ${resolution} 2>&1>${log_path}${tissue}_binomail_diff_analysis_${resolution}.log
    done
done

### TAD 
resolutions=(10000 20000)
for tissue in ${unique_tissue_array[@]}
do
    for resolution in ${resolutions[@]}
    do
        nohup bash ${code_path}TAD/insulation_score.sh ${tissue} ${resolution} 2>&1>${log_path}${tissue}_insulation_score_${resolution}.log &
    done
done

### Loop
resolutions=(10000)
for tissue in ${unique_tissue_array[@]}
do
    for resolution in ${resolutions[@]}
    do
        nohup bash ${code_path}loop/HiCCUP_call_loop.sh ${tissue} ${resolution} 2>&1>${log_path}${tissue}_HiCCUP_call_loop_${resolution}.log &
    done
done

resolutions=(10kb)
for tissue in ${unique_tissue_array[@]}
do
    for resolution in ${resolutions[@]}
    do
        nohup bash ${code_path}loop/mustache_call_loop.sh ${tissue} ${resolution} 2>&1>${log_path}${tissue}_mustache_call_loop_${resolution}.log &
    done
done