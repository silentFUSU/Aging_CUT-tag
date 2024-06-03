tissues=(brain liver testis colon kidney lung spleen muscle pancreas Hip cecum bonemarrow ileum heart thymus stomach skin)
antibodys=(H3K27ac H3K4me1 H3K4me3 H3K27me3 H3K9me3 H3K36me3)
data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/ChromHMM/15_until_skin/state_transfer/E9_to_E1/
mkdir ${result_path}/matrix
mkdir ${result_path}/plot
mkdir ${result_path}/plot/combined_plot
region_lengths=(2000 5000 10000)
for region_length in ${region_lengths[@]}
do
    for tissue in ${tissues[@]}
    do
        for antibody in ${antibodys[@]}
        do
            bw=$(ls ${data_path}${tissue}/${antibody}/bw/*nodup.bw | sort | head -n 2)
            computeMatrix scale-regions -S ${bw} \
                -R ${result_path}/bed/${tissue}.bed \
                --beforeRegionStartLength ${region_length} --startLabel start --endLabel end \
                --regionBodyLength 1000 \
                --afterRegionStartLength ${region_length} \
                --numberOfProcessors 20 \
                --skipZeros -o ${result_path}/matrix/${tissue}_${antibody}_${region_length}_rep1.mat.gz 

            plotProfile -m  ${result_path}/matrix/${tissue}_${antibody}_${region_length}_rep1.mat.gz \
                -out ${result_path}/plot/${tissue}_${antibody}_${region_length}_rep1.pdf --startLabel start --endLabel end \
                --numPlotsPerRow 5 --clusterUsingSamples 1 2 --perGroup --legendLocation best  \
                --plotTitle "${tissue}_${antibody}_${region_length}_rep1"
            
            bw=$(ls ${data_path}${tissue}/${antibody}/bw/*nodup.bw | sort | head -n 4 | tail -n 2)
            computeMatrix scale-regions -S ${bw} \
                -R ${result_path}/bed/${tissue}.bed \
                --beforeRegionStartLength ${region_length} --startLabel start --endLabel end \
                --regionBodyLength 1000 \
                --afterRegionStartLength ${region_length} \
                --numberOfProcessors 20 \
                --skipZeros -o ${result_path}/matrix/${tissue}_${antibody}_${region_length}_rep2.mat.gz 

            plotProfile -m  ${result_path}/matrix/${tissue}_${antibody}_${region_length}_rep2.mat.gz \
                -out ${result_path}/plot/${tissue}_${antibody}_${region_length}_rep2.pdf --startLabel start --endLabel end \
                --numPlotsPerRow 5 --clusterUsingSamples 1 2 --perGroup --legendLocation best  \
                --plotTitle "${tissue}_${antibody}_${region_length}_rep2"
            pdfunite ${result_path}/plot/${tissue}_${antibody}_${region_length}_rep1.pdf ${result_path}/plot/${tissue}_${antibody}_${region_length}_rep2.pdf ${result_path}/plot/combined_plot/${tissue}_${antibody}_${region_length}.pdf

        done
        /storage/zhangyanxiaoLab/suzhuojie/software/pdfjam-3.11/bin/pdfjam ${result_path}/plot/${tissue}_{H3K27ac,H3K4me1,H3K4me3,H3K27me3,H3K9me3,H3K36me3}_${region_length}*.pdf --nup 4x3  --landscape --outfile ${result_path}/plot/combined_plot/${tissue}_${region_length}.pdf
    done
    /storage/zhangyanxiaoLab/suzhuojie/software/pdfjam-3.11/bin/pdfjam ${result_path}/plot/*_{H3K4me3}_${region_length}*.pdf --nup 6x6  --landscape --outfile ${result_path}/plot/combined_plot/all_tissues_${region_length}.pdf
done
