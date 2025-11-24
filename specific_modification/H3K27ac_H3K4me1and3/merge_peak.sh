data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=brain
peaks=$(ls ${data_path}${tissue}/{H3K27ac,H3K4me,H3K4me3}/peaks/macs_narrowpeak/*.narrowPeak)
mkdir ${data_path}${tissue}/H3K27ac_H3K4me1and3
cat $peaks | cut -f 1-3 |sort -k1,1 -k2,2n -|bedtools merge -i stdin | awk -v OFS='\t' 'BEGIN{num=0}{num++; print $0, "peak"num}' - |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}${tissue}/H3K27ac_H3K4me1and3/H3K27ac_H3K4me1and3.bed
bash /storage/zhangyanxiaoLab/suzhuojie/projects/DIPG/code/peaks2matrix/bed_to_saf.sh ${data_path}${tissue}/H3K27ac_H3K4me1and3/H3K27ac_H3K4me1and3.bed ${data_path}${tissue}/H3K27ac_H3K4me1and3/H3K27ac_H3K4me1and3.saf
files=$(ls ${data_path}${tissue}/{H3K27ac,H3K4me,H3K4me3}/bam/*.nodup.bam)
featureCounts -p -a ${data_path}${tissue}/H3K27ac_H3K4me1and3/H3K27ac_H3K4me1and3.saf -o ${data_path}${tissue}/H3K27ac_H3K4me1and3/H3K27ac_H3K4me1and3.counts ${files} -F SAF -T 8 
