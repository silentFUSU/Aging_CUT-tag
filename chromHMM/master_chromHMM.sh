bash chromHMM_bam2bed.sh
bash chromHMM_cellmarkfiletable.sh
bash chromHMM_binarizebed.sh
tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
for tissue in ${tissues[@]}
do
    if [[ "$tissue" == "ovary" || "$tissue" == "mammarygland" || "$tissue" == "uterus" ]]; then
        chrs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX)
    else
        chrs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
    fi
    for chr in ${chrs[@]}
    do
        ln -s /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/combined_analysis_enhancer/binarizedData/${tissue}*${chr}_binary.txt /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/combined_analysis_enhancer/binarizedData_normal_chr/
    done
done

tissues=(aorta BAT bladder bonemarrow brain CB cecum colon heart Hip ileum jejunum kidney liver lung muscle ovary pancreas skin spleen stomach testis thymus tongue uterus mammarygland iWAT)
for tissue in ${tissues[@]}
do
    if [[ "$tissue" == "ovary" || "$tissue" == "mammarygland" || "$tissue" == "uterus" ]]; then
        chrs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX)
    else
        chrs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
    fi
    for chr in ${chrs[@]}
    do
        ln -s /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/combined_analysis_enhancer/binarizedData/${tissue}_*1_${chr}_binary.txt /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/combined_analysis_enhancer/binarizedData_normal_chr_split_rep/rep1/
        ln -s /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/combined_analysis_enhancer/binarizedData/${tissue}_*2_${chr}_binary.txt /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/combined_analysis_enhancer/binarizedData_normal_chr_split_rep/rep2/
    done
done

bash chromHMM_learnmodel.sh
bash chromHMM_comparemodels.sh