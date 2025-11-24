tissues=(colon cecum jejunum ileum)
CUTTAG_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
mkdir ${CUTTAG_path}intestine
antibody=H3K4me1
CUTTAG_files=$(ls ${CUTTAG_path}{colon,cecum,jejunum,ileum}/$antibody/bam/*.bam)
featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_1kb_bins.saf -o ${CUTTAG_path}intestine/${antibody}_1kb_bins.counts ${CUTTAG_files[@]} -F SAF -T 8 

CUTTAG_files=$(ls ${CUTTAG_path}{colon,cecum}/$antibody/bam/*.bam)
featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_1kb_bins.saf -o ${CUTTAG_path}intestine/${antibody}_cecum_colon_1kb_bins.counts ${CUTTAG_files[@]} -F SAF -T 8 

antibodys=(H3K9me3 H3K27me3 H3K36me3)
for antibody in ${antibodys[@]}
do
    CUTTAG_files=$(ls ${CUTTAG_path}{colon,cecum,jejunum,ileum}/$antibody/bam/*.bam)
    featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_10kb_bins.saf -o ${CUTTAG_path}intestine/${antibody}_10kb_bins.counts ${CUTTAG_files[@]} -F SAF -T 8 
done
antibodys=(H3K27ac H3K4me3 H3K4me1)
for antibody in ${antibodys[@]}
do
    CUTTAG_files=$(ls ${CUTTAG_path}{colon,cecum,jejunum,ileum}/$antibody/bam/*.bam)
    featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_1kb_bins.saf -o ${CUTTAG_path}intestine/${antibody}_1kb_bins.counts ${CUTTAG_files[@]} -F SAF -T 8 
done


antibodys=(H3K9me3 H3K27me3 H3K36me3)
for antibody in ${antibodys[@]}
do
    CUTTAG_files=$(ls ${CUTTAG_path}{colon,cecum}/$antibody/bam/*.bam)
    featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_10kb_bins.saf -o ${CUTTAG_path}intestine/${antibody}_cecum_colon_10kb_bins.counts ${CUTTAG_files[@]} -F SAF -T 8 
done

antibodys=(H3K27ac H3K4me3 H3K4me1)
for antibody in ${antibodys[@]}
do
    CUTTAG_files=$(ls ${CUTTAG_path}{colon,cecum}/$antibody/bam/*.bam)
    featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/mm10_1kb_bins.saf -o ${CUTTAG_path}intestine/${antibody}_cecum_colon_1kb_bins.counts ${CUTTAG_files[@]} -F SAF -T 8 
done
