data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=brain
antibody=H3K36me3
ref=mm10
bedtools intersect -a  ${data_path}${tissue}/${antibody}/bed/${antibody}_merge-W1000-G3000-E100.bed -b  ${data_path}${tissue}/${antibody}/bed/${antibody}_merge-W5000-G10000-E100.bed -wa > ${data_path}${tissue}/${antibody}/bed/W1000-G3000-W5000-G10000-intersect.bed

files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}${tissue}/${antibody}/bed/W1000-G3000-W5000-G10000-intersect.bed  ${data_path}${tissue}/${antibody}/bed/W1000-G3000-W5000-G10000-intersect.saf
featureCounts -p -a ${data_path}${tissue}/${antibody}/bed/W1000-G3000-W5000-G10000-intersect.saf -o ${data_path}${tissue}/${antibody}/W1000-G3000-W5000-G10000-intersect.counts ${files} -F SAF -T 8 
