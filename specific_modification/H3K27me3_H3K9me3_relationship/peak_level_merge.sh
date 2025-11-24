data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
antibody1=H3K27me3
antibody2=H3K9me3
tissue=$1
window_size=1000
gap_size=3000
e_value=100
bed1=${data_path}${tissue}/${antibody1}/bed/${antibody1}_young_merge-W${window_size}-G${gap_size}-E${e_value}.bed
bed2=${data_path}${tissue}/${antibody2}/bed/${antibody2}_young_merge-W${window_size}-G${gap_size}-E${e_value}.bed
mkdir ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/
bedtools intersect -a ${bed2} -b ${bed1} -f 0.8   > ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/temp1.bed
bedtools intersect -a ${bed2} -b ${bed1} -F 0.8   > ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/temp2.bed
cat ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/temp1.bed ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/temp2.bed > ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/combined.bed
bedtools sort -i ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/combined.bed | bedtools merge -i - >  ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/${antibody1}_08_${antibody2}_08_merge.bed

bedtools intersect -a ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/${antibody1}_08_${antibody2}_08_merge.bed -b ${bed1} -wao > ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/${antibody1}_intersect_${antibody1}_08_${antibody2}_08_merge.bed

awk '{
    if ($9 != ".") {  # 检查是否有重叠
        lenB = $6 - $5;  # B 文件中对应 peak 的长度
        overlap = $9;  # 交集的长度
        if (lenB != 0) {  # 防止除以零
            percentage = (overlap / lenB) * 100;  # 计算重叠占比
            print $1, $2, $3, $4, $5, $6, percentage;  # 打印结果
        }
    }
}' OFS="\t" ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/${antibody1}_intersect_${antibody1}_08_${antibody2}_08_merge.bed > ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/${antibody1}_percentage_intersect_${antibody1}_08_${antibody2}_08_merge.bed

bedtools intersect -a ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/${antibody1}_08_${antibody2}_08_merge.bed -b ${bed2} -wao > ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/${antibody2}_intersect_${antibody1}_08_${antibody2}_08_merge.bed

awk '{
    if ($9 != ".") {  # 检查是否有重叠
        lenB = $6 - $5;  # B 文件中对应 peak 的长度
        overlap = $9;  # 交集的长度
        if (lenB != 0) {  # 防止除以零
            percentage = (overlap / lenB) * 100;  # 计算重叠占比
            print $1, $2, $3, $4, $5, $6, percentage;  # 打印结果
        }
    }
}' OFS="\t" ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/${antibody2}_intersect_${antibody1}_08_${antibody2}_08_merge.bed > ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/${antibody2}_percentage_intersect_${antibody1}_08_${antibody2}_08_merge.bed

/usr/local/lib64/R/bin/Rscript ~/projects/Aging_CUT_Tag/code/H3K27me3_H3K9me3_relationship/peak_level_intersect_distribution.R ${tissue}

bedtools intersect -a ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/colocalized_tmp.bed -b ${bed2} -wa -wb  | \
awk '{print $1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6}' | sort -k1,1 -k2,2n | bedtools merge > ${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/colocalized.bed

antibodys=(H3K27me3 H3K9me3)
result_path=${data_path}${tissue}/${antibody1}_${antibody2}_intersect/peak_level/
mkdir ${result_path}matrix
mkdir ${result_path}plot

conditions=(colocalized embedded overlap)
for condition in ${conditions[@]}
do
    for antibody in ${antibodys[@]}
    do
        bw=$(ls ${data_path}${tissue}/${antibody}/bw/*nodup.bw | sort | head -n 2)
            computeMatrix scale-regions -S ${bw} \
                -R ${result_path}${condition}.bed \
                --beforeRegionStartLength 2000 --startLabel start --endLabel end \
                --regionBodyLength 2000 \
                --afterRegionStartLength 2000 \
                --numberOfProcessors 20 \
                --skipZeros -o ${result_path}/matrix/${antibody}_${condition}_rep1.mat.gz 

            plotProfile -m  ${result_path}/matrix/${antibody}_${condition}_rep1.mat.gz  \
                -out ${result_path}/plot/${antibody}_${condition}_rep1.pdf --startLabel start --endLabel end \
                --numPlotsPerRow 5 --clusterUsingSamples 1 2 --perGroup --legendLocation best  \
                --plotTitle "${tissue}_${antibody}_${condition}_rep1"
                
        bw=$(ls ${data_path}${tissue}/${antibody}/bw/*nodup.bw | sort | head -n 4 | tail -n 2)
            computeMatrix scale-regions -S ${bw} \
                -R ${result_path}${condition}.bed \
                --beforeRegionStartLength 2000 --startLabel start --endLabel end \
                --regionBodyLength 2000 \
                --afterRegionStartLength 2000 \
                --numberOfProcessors 20 \
                --skipZeros -o ${result_path}/matrix/${antibody}_${condition}_rep2.mat.gz 

            plotProfile -m  ${result_path}/matrix/${antibody}_${condition}_rep2.mat.gz  \
                -out ${result_path}/plot/${antibody}_${condition}_rep2.pdf --startLabel start --endLabel end \
                --numPlotsPerRow 5 --clusterUsingSamples 1 2 --perGroup --legendLocation best  \
                --plotTitle "${tissue}_${antibody}_${condition}_rep2"
    done
    /storage/zhangyanxiaoLab/suzhuojie/software/pdfjam-3.11/bin/pdfjam ${result_path}/plot/*${condition}*.pdf --nup 2x2 --landscape --outfile ${result_path}/plot/${condition}.pdf
done

