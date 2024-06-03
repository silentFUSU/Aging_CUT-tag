data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
# tissues=(brain liver kidney colon testis)
# tissues=(spleen lung muscle pancreas)
tissues=(skin)
ref=mm10
mkdir ${data_path}all/combined_analysis_enhancer/
mkdir ${data_path}all/combined_analysis_enhancer/chromHMMbed
i=1
for tissue in ${tissues[@]}
do
    # mkdir ${data_path}${tissue}/combined_analysis_enhancer/
    # mkdir ${data_path}${tissue}/combined_analysis_enhancer/chromHMMbed

    modifications=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
    for modification in ${modifications[@]}
    do
        samples=$(find ${data_path}${tissue}/${modification}/bam/ -name "*.nodup.bam" -exec basename {} \; | sed 's/\..*//' | sort)
        samples_id=()
        for sample in $samples
        do
            samples_id+=($sample)
        done
        for sample in ${samples_id[@]}
        do
            bedtools bamtobed -i ${data_path}${tissue}/${modification}/bam/${sample}*.nodup.bam > ${data_path}all/combined_analysis_enhancer/chromHMMbed/${sample}.bed &
            # ln -s ${data_path}${tissue}/${modification}/bam/${sample}*.nodup.bam ${data_path}${tissue}/combined_analysis_enhancer/chromHMMbed/${sample}.bam
        done
        j=0
        if [ $i -eq 0 ]; then
            echo -e "${tissue}_young1\t${modification}\t${samples_id[${j}]}.bed" >  ${data_path}all/combined_analysis_enhancer/cellmarkfiletable.txt
            j=$(($j+1))
            echo -e "${tissue}_old1\t${modification}\t${samples_id[${j}]}.bed" >>  ${data_path}all/combined_analysis_enhancer/cellmarkfiletable.txt    
            j=$(($j+1))
            echo -e "${tissue}_young2\t${modification}\t${samples_id[${j}]}.bed" >>  ${data_path}all/combined_analysis_enhancer/cellmarkfiletable.txt 
            j=$(($j+1))
            echo -e "${tissue}_old2\t${modification}\t${samples_id[${j}]}.bed" >> ${data_path}all/combined_analysis_enhancer/cellmarkfiletable.txt
            j=$(($j+1))
            i=1
        else
            echo -e "${tissue}_young1\t${modification}\t${samples_id[${j}]}.bed" >>  ${data_path}all/combined_analysis_enhancer/cellmarkfiletable.txt
            j=$(($j+1))
            echo -e "${tissue}_old1\t${modification}\t${samples_id[${j}]}.bed" >>  ${data_path}all/combined_analysis_enhancer/cellmarkfiletable.txt 
            j=$(($j+1)) 
            echo -e "${tissue}_young2\t${modification}\t${samples_id[${j}]}.bed" >>  ${data_path}all/combined_analysis_enhancer/cellmarkfiletable.txt  
            j=$(($j+1))
            echo -e "${tissue}_old2\t${modification}\t${samples_id[${j}]}.bed" >> ${data_path}all/combined_analysis_enhancer/cellmarkfiletable.txt
            j=$(($j+1))
        fi
    done
    wait
    ChromHMM_path=~/software/ChromHMM/
    # mkdir ${data_path}${tissue}/combined_analysis_enhancer/binarizedData
    java -Xmx32G -jar ${ChromHMM_path}ChromHMM.jar BinarizeBed -b 1000 ${ChromHMM_path}CHROMSIZES/${ref}.txt ${data_path}all/combined_analysis_enhancer/chromHMMbed/ ${data_path}all/combined_analysis_enhancer/cellmarkfiletable.txt ${data_path}all/combined_analysis_enhancer/binarizedData

    state=15
    mkdir ${result_path}${tissue}/ChromHMM/${state}_until_skin
    java -Xmx32G -jar ${ChromHMM_path}ChromHMM.jar LearnModel -b 1000 -p 16 -holdcolumnorder ${data_path}all/combined_analysis_enhancer/binarizedData ${result_path}all/ChromHMM/${state}_until_skin ${state} ${ref}

done