# bash macs2_mergebam_callpeaks_plot.sh -i /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/ -o /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/ -t spleen
func() {
    echo "Usage:"
    echo -e "i:data input path\no:deeptools plot output path\nt:tissue"
}
while getopts ":h:i:o:t:" OPT
do
    case $OPT in
        i) data_path=${OPTARG};;
        o) result_path=${OPTARG};;
        t) tissue=${OPTARG};;        
        h) func;;
        ?) func;;
    esac
done

antibodys=(H3K27ac H3K4me3 H3K4me1)
# antibodys=(H3K27ac H3K4me1)
ref=mm
for antibody in ${antibodys[@]}
do 
    bam=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
    samtools merge -o ${data_path}${tissue}/${antibody}/tmp.merge.bam ${bam} -@ 16
    samtools index ${data_path}${tissue}/${antibody}/tmp.merge.bam -@ 16
    mkdir ${data_path}${tissue}/${antibody}/peaks/
    mkdir ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak

    macs2 callpeak -t ${data_path}${tissue}/${antibody}/tmp.merge.bam  -f BAMPE -n ${antibody} --outdir ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak -g ${ref} --nomodel -q 0.0001  --keep-dup all
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' ${data_path}${tissue}/${antibody}/peaks/macs_narrowpeak/${antibody}_peaks.narrowPeak |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_narrowpeak.bed  
    bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_narrowpeak.bed ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_narrowpeak.saf

    files=$(ls ${data_path}${tissue}/${antibody}/bam/*.nodup.bam)
    featureCounts -p -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_narrowpeak.saf -o ${data_path}${tissue}/${antibody}/${antibody}_macs_narrowpeak.counts ${files} -F SAF -T 8 

    bw=$(ls ${data_path}${tissue}/${antibody}/bw/*nodup.bw)
    mkdir ${data_path}${tissue}/${antibody}/matrix
    computeMatrix scale-regions -S ${bw} \
    -R ${data_path}${tissue}/${antibody}/bed/${antibody}_macs_narrowpeak.bed \
    --beforeRegionStartLength 2000 --startLabel start --endLabel end \
    --regionBodyLength 2000 \
    --afterRegionStartLength 2000 \
    --numberOfProcessors 20 \
    --skipZeros -o ${data_path}${tissue}/${antibody}/matrix/${antibody}_merge_macs2_narrowpeak.mat.gz 
    mkdir ${result_path}${tissue}
    plotProfile -m  ${data_path}${tissue}/${antibody}/matrix/${antibody}_merge_macs2_narrowpeak.mat.gz    \
                -out ${result_path}${tissue}/deeptools/${antibody}_merge_macs2_narrowpeak_each_sample.pdf --startLabel start --endLabel end \
                --numPlotsPerRow 5 --perGroup --legendLocation best \
                --plotTitle "${antibody}"

    plotProfile -m  ${data_path}${tissue}/${antibody}/matrix/${antibody}_merge_macs2_narrowpeak.mat.gz    \
                -out ${result_path}${tissue}/deeptools/${antibody}_merge_macs2_narrowpeak.pdf --startLabel start --endLabel end \
                --numPlotsPerRow 5  --legendLocation best \
                --plotTitle "${antibody}"
                
    rm ${data_path}${tissue}/${antibody}/tmp*
    echo done
done
# data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
# result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
# tissue=brain
# antibody=H3K4me
# ref=mm
# samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
# mkdir ${data_path}${tissue}/${antibody}/peaks
# for sample in ${samples[@]}
# do
# macs2 callpeak -t ${data_path}${tissue}/${antibody}/bam/${sample}.nodup.bam  -f BAMPE -n ${sample} --outdir ${data_path}${tissue}/${antibody}/peaks -g ${ref}  --keep-dup all
# done

# peaks=$(ls ${data_path}${tissue}/${antibody}/peaks/*.narrowPeak)
# mkdir ${data_path}${tissue}/${antibody}/bed/
# cat $peaks | cut -f 1-3 |sort -k1,1 -k2,2n -|bedtools merge -i stdin | awk -v OFS='\t' 'BEGIN{num=0}{num++; print $0, "peak"num}' - |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/npc_yulab/code/bulk_atac/keep_regular_chroms.r > ${data_path}${tissue}/${antibody}/bed/${antibody}_macs2_mergecounts.bed
# bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}${tissue}/${antibody}/bed/${antibody}_macs2_mergecounts.bed ${data_path}${tissue}/${antibody}/bed/${antibody}_macs2_mergecounts.saf
# files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
# featureCounts -p -a ${data_path}${tissue}/${antibody}/bed/${antibody}_macs2_mergecounts.saf -o ${data_path}${tissue}/${antibody}/${antibody}_macs2_mergecounts.counts ${files} -F SAF -T 8 

# bw=$(ls ${data_path}${tissue}/${antibody}/bw/*nodup.bw)
# mkdir ${data_path}${tissue}/${antibody}/matrix
# computeMatrix scale-regions -S ${bw} \
#   -R ${data_path}${tissue}/${antibody}/bed/${antibody}_macs2_mergecounts.bed \
#   --beforeRegionStartLength 2000 --startLabel start --endLabel end \
#   --regionBodyLength 2000 \
#   --afterRegionStartLength 2000 \
#   --numberOfProcessors 20 \
#   --skipZeros -o ${data_path}${tissue}/${antibody}/matrix/${antibody}_macs2.mat.gz 

# plotProfile -m  ${data_path}${tissue}/${antibody}/matrix/${antibody}_macs2.mat.gz   \
#               -out ${result_path}${tissue}/deeptools/${antibody}_macs2.pdf --startLabel start --endLabel end \
#               --numPlotsPerRow 5 --perGroup --legendLocation best \
#               --plotTitle "${antibody}"