data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/MEF_OE_RNA/
# conditions=(MEF_Bmi1 MEF_Cbx2 MEF_Cbx7)
conditions=(MEF_mEzh2 MEF_hEzh2 MEF_mCbx8)
ref=mm10


for condition in ${conditions[@]}
    do
        search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/RNA_search_table.csv
        cleaned_file=$(mktemp)  
        cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
        samples=$(awk -F',' -v c="MEF_Vector" 'NR > 1 && ($1 == c) {print $4}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a sample_array <<<"$samples"  
        echo ${sample_array[@]}
        files=()
        for sample in ${sample_array[@]}
        do
            file=$(ls ${data_path}bam/${sample}*.bam)
            files+=("$file")
        done

        samples=$(awk -F',' -v c="$condition" 'NR > 1 && ($1 == c) {print $4}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a sample_array <<<"$samples"  
        echo ${sample_array[@]}
        for sample in ${sample_array[@]}
        do
            file=$(ls ${data_path}bam/${sample}*.bam)
            files+=("$file")
        done
        featureCounts -a /storage/zhangyanxiaoLab/share/gtf/mm10.gencode.vM25.annotation.gtf -p -o ${data_path}MEF_Vector_${condition}_merge.counts ${files[@]} -F GTF -T 10 -t exon -g gene_name 
    done

files=$(ls ${data_path}bam/*.bam)
featureCounts -a /storage/zhangyanxiaoLab/share/gtf/mm10.gencode.vM25.annotation.gtf -p -o ${data_path}all_merge.counts ${files[@]} -F GTF -T 10 -t exon -g gene_name 