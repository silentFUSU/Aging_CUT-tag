data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF_OE/
antibodys=(H3K27me3 H2AK119ub1)
conditions=(MEF_Vector MEF_Bmi1 MEF_Cbx2 MEF_Cbx7)
ref=mm10
for antibody in ${antibodys[@]}
do
    for condition in ${conditions[@]}
    do
        search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff_batch.csv
        cleaned_file=$(mktemp)  
        cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
        samples=$(awk -F',' -v c="$condition" -v a="$antibody" 'NR > 1 && ($1 == c) && ($2 == a) {print $3}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a sample_array <<<"$samples"  
        echo ${sample_array[@]}
        files=()
        for sample in ${sample_array[@]}
        do
            file=$(ls ${data_path}${antibody}/${condition}/bam/${sample}.bam)
            files+=("$file")
        done
        echo "samtools merge -o ${data_path}${antibody}/${condition}/tmp.merge.bam ${files[@]} -@ 16"
        samtools merge -f -o ${data_path}${antibody}/${condition}/tmp.merge.bam ${files[@]} -@ 16
        samtools index ${data_path}${antibody}/${condition}/tmp.merge.bam -@ 16
        window_size=5000
        gap_size=10000
        e_value=100
        mkdir -p ${data_path}${antibody}/${condition}/peaks/
        sicer  -t ${data_path}${antibody}/${condition}/tmp.merge.bam -o ${data_path}${antibody}/${condition}/peaks/ -s ${ref} -w ${window_size} -rt 16 -f 300 -egf 0.8 -fdr 0.01 -g ${gap_size} -e ${e_value} -cpu 21
        rm ${data_path}${antibody}/${condition}/tmp*
        mkdir -p ${data_path}${antibody}/${condition}/bed/
        awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR, $4}' ${data_path}${antibody}/${condition}/peaks/tmp.merge-W${window_size}-G${gap_size}.scoreisland |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r >  ${data_path}${antibody}/${condition}/bed/${antibody}_${condition}_merge-W${window_size}-G${gap_size}-E${e_value}.bed  
    
    done
done

conditions=(MEF_Bmi1 MEF_Cbx2 MEF_Cbx7)
for antibody in ${antibodys[@]}
do
    for condition in ${conditions[@]}
    do
    cat  ${data_path}${antibody}/${condition}/bed/${antibody}_MEF_Vector_merge-W${window_size}-G${gap_size}-E${e_value}.bed   ${data_path}${antibody}/${condition}/bed/${antibody}_${condition}_merge-W${window_size}-G${gap_size}-E${e_value}.bed | \
        sort -k1,1 -k2,2n | \
        bedtools merge  >  ${data_path}${antibody}/${antibody}_MEF_Vector_${condition}_merge-W${window_size}-G${gap_size}-E${e_value}.bed
    done
done