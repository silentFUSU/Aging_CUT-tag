tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
antibodys=(H3K27me3 H3K9me3 H3K36me3)
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
        search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff_batch.csv
        cleaned_file=$(mktemp)  
        cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
        samples=$(awk -F',' -v t="$tissue" -v y="3m" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == y) && ($2 == a) {print $3}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a young_array <<<"$samples"  
        echo ${young_array[@]}
        young=()
        for sample in ${young_array[@]}
        do
        file=$(ls ${data_path}${tissue}/${antibody}/bam/${sample}*.bam)
        young+=("$file")
        done
        echo ${young[@]}
        samples=$(awk -F',' -v t="$tissue" -v y="24m" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == y) && ($2 == a) {print $3}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a old_array <<<"$samples"
        echo ${old_array[@]}
        old=()
        for sample in ${old_array[@]}
        do
        file=$(ls ${data_path}${tissue}/${antibody}/bam/${sample}*.bam)
        old+=("$file")
        done
        echo ${old[@]}
        echo "samtools merge -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young[@]} -@ 16"
        echo "samtools merge -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old[@]} -@ 16"
        samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young[@]} -@ 16 &
        samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old[@]} -@ 16 &
        wait
        samtools index ${data_path}${tissue}/${antibody}/tmp.young.merge.bam -@ 16 &
        samtools index ${data_path}${tissue}/${antibody}/tmp.old.merge.bam -@ 16 &
        wait
        # binsize=1000
        # smoothLength=3000
        binsize=10000
        smoothLength=100000
        bamCoverage -p 30 -e 100 --binSize ${binsize} --smoothLength ${smoothLength} -b ${data_path}${tissue}/${antibody}/tmp.young.merge.bam -o ${data_path}${tissue}/${antibody}/bw/young_bs${binsize}.bw --normalizeUsing RPKM &
        bamCoverage -p 30 -e 100 --binSize ${binsize} --smoothLength ${smoothLength} -b ${data_path}${tissue}/${antibody}/tmp.old.merge.bam -o ${data_path}${tissue}/${antibody}/bw/old_bs${binsize}.bw --normalizeUsing RPKM &
        wait
        rm ${data_path}${tissue}/${antibody}/tmp*
    done
done

antibodys=(H3K27ac H3K4me1 H3K4me3)
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
        search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff_batch.csv
        cleaned_file=$(mktemp)  
        cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
        samples=$(awk -F',' -v t="$tissue" -v y="3m" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == y) && ($2 == a) {print $3}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a young_array <<<"$samples"  
        echo ${young_array[@]}
        young=()
        for sample in ${young_array[@]}
        do
        file=$(ls ${data_path}${tissue}/${antibody}/bam/${sample}*.bam)
        young+=("$file")
        done
        echo ${young[@]}
        samples=$(awk -F',' -v t="$tissue" -v y="24m" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == y) && ($2 == a) {print $3}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a old_array <<<"$samples"
        echo ${old_array[@]}
        old=()
        for sample in ${old_array[@]}
        do
        file=$(ls ${data_path}${tissue}/${antibody}/bam/${sample}*.bam)
        old+=("$file")
        done
        echo ${old[@]}
        echo "samtools merge -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young[@]} -@ 16"
        echo "samtools merge -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old[@]} -@ 16"
        samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young[@]} -@ 16 &
        samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old[@]} -@ 16 &
        wait
        samtools index ${data_path}${tissue}/${antibody}/tmp.young.merge.bam -@ 16 &
        samtools index ${data_path}${tissue}/${antibody}/tmp.old.merge.bam -@ 16 &
        wait
        binsize=50
        bamCoverage -p 30 -e 100 --binSize ${binsize} -b ${data_path}${tissue}/${antibody}/tmp.young.merge.bam -o ${data_path}${tissue}/${antibody}/bw/young_bs${binsize}.bw --normalizeUsing RPKM &
        bamCoverage -p 30 -e 100 --binSize ${binsize} -b ${data_path}${tissue}/${antibody}/tmp.old.merge.bam -o ${data_path}${tissue}/${antibody}/bw/old_bs${binsize}.bw --normalizeUsing RPKM &
        wait
        rm ${data_path}${tissue}/${antibody}/tmp*
    done
done


antibodys=(ATAC)
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ATAC/
        search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/ATAC_search_table_batch.csv
        cleaned_file=$(mktemp)  
        cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
        samples=$(awk -F',' -v t="$tissue" -v y="3m" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == y) && ($2 == a) {print $3}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a young_array <<<"$samples"  
        echo ${young_array[@]}
        young=()
        for sample in ${young_array[@]}
        do
        file=$(ls ${data_path}${tissue}/${antibody}/bam/${sample}*.bam)
        young+=("$file")
        done
        echo ${young[@]}
        samples=$(awk -F',' -v t="$tissue" -v y="24m" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == y) && ($2 == a) {print $3}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a old_array <<<"$samples"
        echo ${old_array[@]}
        old=()
        for sample in ${old_array[@]}
        do
        file=$(ls ${data_path}${tissue}/${antibody}/bam/${sample}*.bam)
        old+=("$file")
        done
        echo ${old[@]}
        echo "samtools merge -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young[@]} -@ 16"
        echo "samtools merge -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old[@]} -@ 16"
        samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.young.merge.bam ${young[@]} -@ 16 &
        samtools merge -f -o ${data_path}${tissue}/${antibody}/tmp.old.merge.bam ${old[@]} -@ 16 &
        wait
        samtools index ${data_path}${tissue}/${antibody}/tmp.young.merge.bam -@ 16 &
        samtools index ${data_path}${tissue}/${antibody}/tmp.old.merge.bam -@ 16 &
        wait
        binsize=50
        bamCoverage -p 30 -e 100 --binSize ${binsize} -b ${data_path}${tissue}/${antibody}/tmp.young.merge.bam -o ${data_path}${tissue}/${antibody}/bw/young_bs${binsize}.bw --normalizeUsing RPKM &
        bamCoverage -p 30 -e 100 --binSize ${binsize} -b ${data_path}${tissue}/${antibody}/tmp.old.merge.bam -o ${data_path}${tissue}/${antibody}/bw/old_bs${binsize}.bw --normalizeUsing RPKM &
        wait
        rm ${data_path}${tissue}/${antibody}/tmp*
    done
done

antibodys=(RNA)
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
        search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/RNA_search_table.csv
        cleaned_file=$(mktemp)  
        cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
        samples=$(awk -F',' -v t="$tissue" -v y="3m" -v a="$antibody" 'NR > 1 && ($2 == t) && ($6 == y) && ($3 == a) {print $4}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a young_array <<<"$samples"  
        echo ${young_array[@]}
        young=()
        for sample in ${young_array[@]}
        do
        file=$(ls ${data_path}${tissue}/bam/${sample}*.bam)
        young+=("$file")
        done
        echo ${young[@]}
        samples=$(awk -F',' -v t="$tissue" -v y="24m" -v a="$antibody" 'NR > 1 && ($2 == t) && ($6 == y) && ($3 == a) {print $4}' "$cleaned_file")
        IFS=$'\n' read -rd '' -a old_array <<<"$samples"
        echo ${old_array[@]}
        old=()
        for sample in ${old_array[@]}
        do
        file=$(ls ${data_path}${tissue}/bam/${sample}*.bam)
        old+=("$file")
        done
        echo ${old[@]}
        echo "samtools merge -o ${data_path}${tissue}/tmp.young.merge.bam ${young[@]} -@ 16"
        echo "samtools merge -o ${data_path}${tissue}/tmp.old.merge.bam ${old[@]} -@ 16"
        samtools merge -f -o ${data_path}${tissue}/tmp.young.merge.bam ${young[@]} -@ 16 &
        samtools merge -f -o ${data_path}${tissue}/tmp.old.merge.bam ${old[@]} -@ 16 &
        wait
        samtools index ${data_path}${tissue}/tmp.young.merge.bam -@ 16 &
        samtools index ${data_path}${tissue}/tmp.old.merge.bam -@ 16 &
        wait
        binsize=50
        bamCoverage -p 30 -e 100 --binSize ${binsize} -b ${data_path}${tissue}/tmp.young.merge.bam -o ${data_path}${tissue}/bw/young_bs${binsize}.bw --normalizeUsing RPKM &
        bamCoverage -p 30 -e 100 --binSize ${binsize} -b ${data_path}${tissue}/tmp.old.merge.bam -o ${data_path}${tissue}/bw/old_bs${binsize}.bw --normalizeUsing RPKM &
        wait
        rm ${data_path}${tissue}/tmp*
    done
done