data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
# tissues=(brain liver kidney colon testis)
# tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum bonemarrow)
tissues=(bladder aorta tongue)
antibodys=(H3K36me3 H3K27me3 H3K9me3)
ref=mm10
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do 
        samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
        binsize=1000
        smoothLength=3000
        for sample in ${samples[@]}
        do
            bamCoverage -p 30 -e 100 --binSize ${binsize} --smoothLength ${smoothLength} -b ${data_path}${tissue}/${antibody}/bam/${sample}.nodup.bam -o ${data_path}${tissue}/${antibody}/bw/${sample}_bs${binsize}.bw --normalizeUsing RPKM
        done
    done
done

for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do 
        samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
        binsize=10000
        smoothLength=100000
        for sample in ${samples[@]}
        do
            bamCoverage -p 30 -e 100 --binSize ${binsize} --smoothLength ${smoothLength} -b ${data_path}${tissue}/${antibody}/bam/${sample}.nodup.bam -o ${data_path}${tissue}/${antibody}/bw/${sample}_bs${binsize}.bw --normalizeUsing RPKM
        done
    done
done

# data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/
# tissues=(cellular_aging_GSE133292)
# antibodys=(Chip_seq)
# ref=hg38
# mkdir ${data_path}${tissue}/${antibody}/bw/
# for tissue in ${tissues[@]}
# do
#     for antibody in ${antibodys[@]}
#     do 
#         samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
#         binsize=10000
#         smoothLength=100000
#         for sample in ${samples[@]}
#         do
#             bamCoverage -p 30 -e 100 --binSize ${binsize} --smoothLength ${smoothLength} -b ${data_path}${tissue}/${antibody}/bam/${sample}.nodup.bam -o ${data_path}${tissue}/${antibody}/bw/${sample}_bs${binsize}.bw --normalizeUsing RPKM
#         done
#     done
# done