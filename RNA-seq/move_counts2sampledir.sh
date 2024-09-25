data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240918_LLX727_RNA/
target_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
samples=$(find ${data_path}featureCounts -type f -name "*.counts" -exec basename {} \; | sed 's/\.counts//' | sort)  
tissues=(ovary)
# tissues=(jejunum)
# tissues=(cecum cecum cecum cecum colon colon colon colon)
i=0
for sample in ${samples[@]} 
do  
    mv ${data_path}featureCounts/${sample}* ${target_path}${tissues[i]}/counts/
    mv ${data_path}bam/${sample}* ${target_path}${tissues[i]}/bam/
    mv ${data_path}bigWig/${sample}* ${target_path}${tissues[i]}/bw/
    ((i++))
done  