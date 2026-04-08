data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF_OE/
antibodys=(H3K27me3 H2AK119ub1)
conditions=(MEF_Vector MEF_Bmi1 MEF_Cbx2 MEF_Cbx7)
ref=mm10
for antibody in ${antibodys[@]}
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
        binsize=5000
        smoothLength=10000
        bamCoverage -p 20 -e 100 --binSize ${binsize} --smoothLength ${smoothLength} -b ${data_path}${antibody}/${condition}/tmp.merge.bam -o ${data_path}${antibody}/${condition}/bw/${antibody}_${condition}.bw --normalizeUsing RPKM &
done
