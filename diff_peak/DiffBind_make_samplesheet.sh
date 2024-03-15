data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=brain
antibody=H3K27me3

echo "SampleID,Factor,Replicate,bamReads,Peaks" > ${data_path}${tissue}/${antibody}/SampleSheet.csv
samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
sample1=$(echo $samples | awk '{print $1}')  
sample2=$(echo $samples | awk '{print $2}')  
sample3=$(echo $samples | awk '{print $3}')  
sample4=$(echo $samples | awk '{print $4}')
echo "young1,young,1,${data_path}${tissue}/${antibody}/bam/${sample1}.nodup.bam,${data_path}${tissue}/${antibody}/peaks/tmp.merge-W1000-G3000.scoreisland" > ${data_path}${tissue}/${antibody}/SampleSheet.csv
echo "old1,old,1,${data_path}${tissue}/${antibody}/bam/${sample2}.nodup.bam,${data_path}${tissue}/${antibody}/peaks/tmp.merge-W1000-G3000.scoreisland" >> ${data_path}${tissue}/${antibody}/SampleSheet.
echo "young2,young,2,${data_path}${tissue}/${antibody}/bam/${sample3}.nodup.bam,${data_path}${tissue}/${antibody}/peaks/tmp.merge-W1000-G3000.scoreisland" >> ${data_path}${tissue}/${antibody}/SampleSheet.csv
echo "old2,old,2,${data_path}${tissue}/${antibody}/bam/${sample4}.nodup.bam,${data_path}${tissue}/${antibody}/peaks/tmp.merge-W1000-G3000.scoreisland" >> ${data_path}${tissue}/${antibody}/SampleSheet.csv