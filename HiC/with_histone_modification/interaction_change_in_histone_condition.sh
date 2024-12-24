data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=$1
ref=mm10
# antibodys=(H3K27me3 H3K9me3 H3K36me3 H3K4me1 H3K4me3 H3K27ac)
antibodys=(H3K9me3 H3K27me3)
for antibody in ${antibodys[@]}
do
    search_table=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/CUTTag_search_table_used_in_diff.csv
    cleaned_file=$(mktemp)  
    cat "$search_table" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  
    young_samples=$(awk -F',' -v t="$tissue" -v y="3m" -v a="$antibody" 'NR > 1 && ($1 == t) && ($5 == y) && ($2 == a) {print $3}' "$cleaned_file")
    IFS=$'\n' read -rd '' -a young_array <<<"$young_samples"  
    mkdir -p ${data_path}HiC/${tissue}/with_histone
    mkdir -p ${data_path}HiC/${tissue}/with_histone/${antibody}
    mkdir -p ${data_path}HiC/${tissue}/with_histone/${antibody}/peaks
    for sample in ${young_samples[@]}
    do
        macs2 callpeak -B --SPMR --broad --nomodel --nolambda -t ${data_path}${tissue}/${antibody}/bam/${sample}*.nodup.bam  -f BAMPE -n ${sample} --outdir ${data_path}HiC/${tissue}/with_histone/${antibody}/peaks -g mm -q 0.0001  --keep-dup all &
    done
    wait

    for sample in ${young_samples[@]}
    do
        python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/HiC/with_histone_modification/bedtobins.py -i ${data_path}HiC/${tissue}/with_histone/${antibody}/peaks/${sample}_treat_pileup.bdg -o ${data_path}HiC/${tissue}/with_histone/${antibody}/${sample}_10kb.txt &
    done
    wait
done