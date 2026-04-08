data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF_OE/
antibodys=(H3K27me3 H2AK119ub1)
conditions=(MEF_Bmi1 MEF_Cbx2 MEF_Cbx7)
window_size=5000
gap_size=10000
e_value=100
for antibody in ${antibodys[@]}
do
    for condition in ${conditions[@]}
    do
        search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff_batch.csv
        cleaned_file=$(mktemp)  
        cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
        samples=$(awk -F',' -v c="MEF_Vector" -v a="$antibody" 'NR > 1 && ($1 == c) && ($2 == a) {print $3}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a sample_array <<<"$samples"  
        echo ${sample_array[@]}
        files=()
        for sample in ${sample_array[@]}
        do
            file=$(ls ${data_path}${antibody}/MEF_Vector/bam/${sample}.bam)
            files+=("$file")
        done

        samples=$(awk -F',' -v c="$condition" -v a="$antibody" 'NR > 1 && ($1 == c) && ($2 == a) {print $3}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a sample_array <<<"$samples"  
        echo ${sample_array[@]}
        for sample in ${sample_array[@]}
        do
            file=$(ls ${data_path}${antibody}/${condition}/bam/${sample}.bam)
            files+=("$file")
        done
        bed=${data_path}${antibody}/${antibody}_MEF_Vector_${condition}_merge-W${window_size}-G${gap_size}-E${e_value}.bed
        saf=${data_path}${antibody}/${antibody}_MEF_Vector_${condition}_merge-W${window_size}-G${gap_size}-E${e_value}.saf
        bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${bed} ${saf}
        featureCounts -p -a ${saf} -o ${data_path}${antibody}/${antibody}_MEF_Vector_${condition}_merge-W${window_size}-G${gap_size}-E${e_value}.counts ${files[@]} -F SAF -T 8 
    done
done