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
        for sample in ${sample_array[@]}
        do
            sf=$(grep $sample /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data_transposon2/20260211_DYQ_CUTTag_SAMPLE/all_sample.qc.txt |cut -f 8)
            echo $sample $sf
            bamCoverage --scaleFactor $sf -b ${data_path}${antibody}/${condition}/bam/${sample}.bam -o ${data_path}${antibody}/${condition}/bw/${sample}.scaled.filt.srt2.bw --outFileFormat bigwig --binSize 1000 --smoothLength 3000 --numberOfProcessors 6
        done
    done
done