data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
# tissues=(brain liver kidney colon testis)
# tissues=(spleen lung muscle pancreas)
tissues=(lung)
for tissue in ${tissues[@]}
do
    ref=mm10
    samples=$(find ${data_path}${tissue}/*/bam/ -name "*.nodup.bam" -exec basename {} \; | sed 's/\..*//')
    samples_id=()
    for sample in $samples
    do
        extracted=$(echo $sample | cut -d '_' -f 1)
        samples_id+=($extracted)
    done
    mkdir ${data_path}${tissue}/combined_analysis_enhancer/
    mkdir ${data_path}${tissue}/combined_analysis_enhancer/chromHMMbed
    for sample in ${samples_id[@]}
    do
        bedtools bamtobed -i ${data_path}${tissue}/*/bam/${sample}*.nodup.bam > ${data_path}${tissue}/combined_analysis_enhancer/chromHMMbed/${sample}.bed &
    done
    wait
    modifications=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
    i=0
    for modification in ${modifications[@]}
    do
        if [ $i -eq 0 ]; then
            echo -e "young1\t${modification}\t${samples_id[${i}]}.bed" >  ${data_path}${tissue}/combined_analysis_enhancer/chromHMMbed/cellmarkfiletable.txt
            i=$((i+1))    
            echo -e "old1\t${modification}\t${samples_id[${i}]}.bed" >>  ${data_path}${tissue}/combined_analysis_enhancer/chromHMMbed/cellmarkfiletable.txt
            i=$((i+1))    
            echo -e "young2\t${modification}\t${samples_id[${i}]}.bed" >>  ${data_path}${tissue}/combined_analysis_enhancer/chromHMMbed/cellmarkfiletable.txt
            i=$((i+1))    
            echo -e "old2\t${modification}\t${samples_id[${i}]}.bed" >> ${data_path}${tissue}/combined_analysis_enhancer/chromHMMbed/cellmarkfiletable.txt
            i=$((i+1))    
        else
            echo -e "young1\t${modification}\t${samples_id[${i}]}.bed" >>  ${data_path}${tissue}/combined_analysis_enhancer/chromHMMbed/cellmarkfiletable.txt
            i=$((i+1))    
            echo -e "old1\t${modification}\t${samples_id[${i}]}.bed" >>  ${data_path}${tissue}/combined_analysis_enhancer/chromHMMbed/cellmarkfiletable.txt
            i=$((i+1))    
            echo -e "young2\t${modification}\t${samples_id[${i}]}.bed" >>  ${data_path}${tissue}/combined_analysis_enhancer/chromHMMbed/cellmarkfiletable.txt
            i=$((i+1))    
            echo -e "old2\t${modification}\t${samples_id[${i}]}.bed" >> ${data_path}${tissue}/combined_analysis_enhancer/chromHMMbed/cellmarkfiletable.txt
            i=$((i+1))  
        fi
    done

    ChromHMM_path=~/software/ChromHMM/
    mkdir ${data_path}${tissue}/combined_analysis_enhancer/binarizedData
    java -Xmx16G -jar ${ChromHMM_path}ChromHMM.jar BinarizeBed -b 1000 ${ChromHMM_path}CHROMSIZES/${ref}.txt ${data_path}${tissue}/combined_analysis_enhancer/chromHMMbed \
        ${data_path}${tissue}/combined_analysis_enhancer/chromHMMbed/cellmarkfiletable.txt ${data_path}${tissue}/combined_analysis_enhancer/binarizedData

    state=15
    mkdir ${result_path}${tissue}/ChromHMM/${state}
    java -Xmx16G -jar ${ChromHMM_path}ChromHMM.jar LearnModel -p 8 -holdcolumnorder ${data_path}${tissue}/combined_analysis_enhancer/binarizedData ${result_path}${tissue}/ChromHMM/${state} ${state} ${ref}

done