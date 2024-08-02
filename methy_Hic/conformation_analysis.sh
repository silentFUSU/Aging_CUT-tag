data_path=~/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/WJH_Mousebrain_C1_BSseq/WJH_Mousebrain_C1_BSseq_101/
## FanC Compartment analysis
fanc compartments ${data_path}aligned/inter_30.hic@1mb ${data_path}architecture/compartments/fanc_1mb.ab 2>&1>${data_path}fanc_compartments.log &

fanc compartments -g /storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.fa \
                  -d ${data_path}architecture/compartments/fanc_1mb.domains_gc.bed \
                  ${data_path}architecture/compartments/fanc_1mb.ab 2>&1>${data_path}fanc_compartments_generate_bed.log &

fanc compartments -v ${data_path}architecture/compartments/fanc_1mb.ev.txt \
                  ${data_path}architecture/compartments/fanc_1mb.ab

fancplot -o ${data_path}architecture/compartments/fanc_1mb.ab.png chr1 \
     -p square ${data_path}architecture/compartments/fanc_1mb.ab \
     -vmin -0.75 -vmax 0.75 -c RdBu_r \
     -p line ${data_path}architecture/compartments/fanc_1mb.ev.txt

## FanC TAD analysis
fancplot -o ${data_path}architecture/TAD/FanC/fanc_100kb_tads.png chr1:4mb-20mb \
     -p triangular ${data_path}aligned/inter_30.hic@100kb  -m 10000000  -vmin 0 -vmax 100

fancplot -o ${data_path}architecture/TAD/FanC/fanc_50kb_tads.png chr1:4mb-20mb \
     -p triangular ${data_path}aligned/inter_30.hic@50kb  -m 10000000  -vmin 0 -vmax 40 

fanc insulation ${data_path}aligned/inter_30.hic@100kb \
     ${data_path}architecture/TAD/fanc_100kb.insulation \
     -w 1000000 1500000 2000000 2500000 3000000 3500000 4000000

fanc insulation ${data_path}architecture/TAD/fanc_100kb.insulation \
                -o bed

fanc insulation ${data_path}architecture/TAD/fanc_100kb.insulation

fancplot -o ${data_path}architecture/TAD/FanC/fanc_100kb_tads_insulation.png  chr1:4mb-20mb \
     -p triangular ${data_path}aligned/inter_30.hic@100kb -m 10000000  -vmin 0 -vmax 100  \
     -p scores ${data_path}architecture/TAD/FanC/fanc_100kb.insulation

fancplot --width 6 -o ${data_path}architecture/TAD/fanc_100kb_tads_insulation_1mb.png \
               chr18:18mb-28mb \
               -p triangular ${data_path}aligned/inter_30.hic@100kb -m 4000000 \
               -vmin 0 -vmax 500 \
               -p line ${data_path}architecture/TAD/fanc_100kb.insulation_1mb.bed \
               ${data_path}architecture/TAD/fanc_100kb.insulation_2mb.bed \
               -l "1mb" "2mb"
## FanC PCA analysis
fanc pca -n "BS_101" "BS_102" "C1_101" "C1_102" \
         -Z -s 1000000 -f -p ~/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/FanC.pca.png \
         ~/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/WJH_Mousebrain_BSseq/WJH_Mousebrain_BSseq_101/aligned/inter_30.hic@1mb \
         ~/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/WJH_Mousebrain_BSseq/WJH_Mousebrain_BSseq_102/aligned/inter_30.hic@1mb \
         ~/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/WJH_Mousebrain_C1_BSseq/WJH_Mousebrain_C1_BSseq_101/aligned/inter_30.hic@1mb \
         ~/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/WJH_Mousebrain_C1_BSseq/WJH_Mousebrain_C1_BSseq_102/aligned/inter_30.hic@1mb  \
         ~/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/lowc.pca

## Homer call compartment

awk 'BEGIN { OFS = "\t" } {  
  strand1 = ($1 == 0 ? "+" : ($1 == 16 ? "-" : $1))  
  strand2 = ($7 == 0 ? "+" : ($7 == 16 ? "-" : $7))  
  print "pairs" NR, $2, strand1, $3, $6, $5, strand2  
}' ${data_path}aligned/merged30.txt > ${data_path}architecture/compartments/Homer/merged30_homer.txt 

# makeTagDirectory ${data_path}architecture/compartments/Homer/HiCsummary/ -format HiCsummary ${data_path}architecture/compartments/Homer/merged30_homer.txt 
tagDir2hicFile.pl ${data_path}architecture/compartments/Homer/HiCsummary/ -juicer ${data_path}aligned/inter_30.hic  -genome mm10 -juicerExe ~/software/juicer/scripts/common/juicer_tools  -p 10
runHiCpca.pl ${data_path}architecture/compartments/Homer/Homer_1m ${data_path}architecture/compartments/Homer/HiCsummary/ -res 1000000 -window 2000000 -cpu 4 -std 1 -corrDepth 1 -min 0 -genome mm10


# outputdir=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/WJH_Mousebrain_C1_BSseq/WJH_Mousebrain_C1_BSseq_101/aligned

# ~/software/juicer/scripts/common/index_by_chr.awk ${outputdir}/merged30.txt 500000 > ${outputdir}/merged30_index.txt

# ~/software/juicer/scripts/common/juicer_tools pre -n -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m \
#      -f /storage/zhangyanxiaoLab/suzhuojie/software/juicer/restriction_sites/mm10_DpnII.txt -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100 --threads 8 -i ${outputdir}/merged30_index.txt \
#      -t ${outputdir}/HIC30_tmp ${outputdir}/merged30.txt ${outputdir}/inter_30.hic /storage/zhangyanxiaoLab/suzhuojie/software/juicer/references/mm10/mm10.chrom.sizes
# nohup ~/software/juicer/scripts/common/juicer_tools addNorm --threads 8 ${outputdir}/inter_30.hic 2>&1>/dev/null &

# hicConvertFormat -m ${data_path}aligned/inter_30.hic -o ${data_path}architecture/compartments/Homer/inter_30.homer --inputFormat hic --outputFormat homer 
data_path2=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240725_methylHiC/lambda_test_all_data/WJH_Mousebrain_C1_BSseq/WJH_Mousebrain_C1_BSseq_101/
mkdir ${data_path2}architecture
mkdir ${data_path2}architecture/compartments
mkdir ${data_path2}architecture/compartments/HiCExplorer/
hic2cool convert ${data_path2}aligned/inter_30.hic ${data_path2}architecture/compartments/HiCExplorer/inter_30.cool -r 0
hic2cool convert ${data_path2}aligned/inter_30.hic ${data_path2}architecture/compartments/HiCExplorer/inter_30_5000.cool -r 5000 
hicConvertFormat -m ${data_path2}architecture/compartments/HiCExplorer/inter_30_5000.cool --inputFormat cool --outputFormat h5 -o ${data_path2}architecture/compartments/HiCExplorer/inter_30_5000.h5 --correction_name KR
hicConvertFormat -m ${data_path2}architecture/compartments/HiCExplorer/inter_30_5000.cool --inputFormat cool --outputFormat homer -o ${data_path2}architecture/compartments/HiCExplorer/inter_30_5000.homer --correction_name KR

nohup java -jar -Xmx48000m  -Djava.awt.headless=true -jar /storage/zhangyanxiaoLab/suzhuojie/software/juicer/scripts/common/juicer_tools.jar arrowhead \
    --threads 8 -k KR -m 2000 -r 5000 ${data_path}aligned/inter_30.hic ${data_path}architecture/TAD/juicer_TAD_5000 2>&1>>${data_path}juicer_TAD.log &