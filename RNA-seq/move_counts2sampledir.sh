data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240524_LLX2/
target_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
samples=$(find ${data_path}featureCounts -type f -name "*.counts" -exec basename {} \; | sed 's/\.counts//' | sort)  
tissues=(jejunum jejunum testis testis testis testis liver liver liver liver)
# tissues=(jejunum)
# tissues=(cecum cecum cecum cecum colon colon colon colon)
i=0
for sample in ${samples[@]} 
do  
    ln -s ${data_path}featureCounts/${sample}* ${target_path}${tissues[i]}/counts/
    ln -s ${data_path}bam/${sample}* ${target_path}${tissues[i]}/bam/
    ln -s ${data_path}bigWig/${sample}* ${target_path}${tissues[i]}/bw/
    ((i++))
done  