data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
antibodys=(H3K9me3)
ref=mm10
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
tissues=$1
for tissue in ${tissues[@]}
do
    for antibody in ${antibodys[@]}
    do
        samples=$(find ${data_path}${tissue}/${antibody}/bam/ -name "*.bam" -exec basename {} \; | sed 's/\..*//')
        files=$(ls ${data_path}${tissue}/${antibody}/bam/*.bam)
        featureCounts -p -a /storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS_100kb.saf -o ${data_path}${tissue}/${antibody}/${antibody}_gene_TSS_100kb.counts ${files} -F SAF -T 8 
    done
done