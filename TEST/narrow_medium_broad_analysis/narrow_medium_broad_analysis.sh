data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
antibody=H3K27me3
bar_low=10000
bar_high=1000000
tissue=liver
window_size=100000
gap_size=300000
# window_size=5000
# gap_size=10000
e_value=100
awk -F',' -v bar_low="${bar_low}" 'NR>1 && $6 <= bar_low {print $2 "\t" $3 "\t" $4}' ${data_path}${tissue}/${antibody}/${antibody}_merge-W${window_size}-G${gap_size}-E${e_value}_diff.csv | tr -d "\"" >  ${data_path}${tissue}/${antibody}/bed/less_than_${bar_low}-W${window_size}-G${gap_size}-E${e_value}.bed  
awk -F',' -v bar_low="${bar_high}" 'NR>1 && $6 >= bar_low {print $2 "\t" $3 "\t" $4}' ${data_path}${tissue}/${antibody}/${antibody}_merge-W${window_size}-G${gap_size}-E${e_value}_diff.csv | tr -d "\"" >  ${data_path}${tissue}/${antibody}/bed/higher_than_${bar_high}-W${window_size}-G${gap_size}-E${e_value}.bed  
bw=$(ls ${data_path}${tissue}/${antibody}/bw/*nodup.bw)
computeMatrix scale-regions -S ${bw} \
    -R ${data_path}${tissue}/${antibody}/bed/less_than_${bar_low}-W${window_size}-G${gap_size}-E${e_value}.bed \
    --beforeRegionStartLength 2000 --startLabel start --endLabel end \
    --regionBodyLength 2000 \
    --afterRegionStartLength 2000 \
    --numberOfProcessors 20 \
    --skipZeros -o ${data_path}${tissue}/${antibody}/matrix/${antibody}_less_than_${bar_low}-W${window_size}-G${gap_size}-E${e_value}.mat.gz &

computeMatrix scale-regions -S ${bw} \
    -R ${data_path}${tissue}/${antibody}/bed/higher_than_${bar_high}-W${window_size}-G${gap_size}-E${e_value}.bed \
    --beforeRegionStartLength 2000 --startLabel start --endLabel end \
    --regionBodyLength 2000 \
    --afterRegionStartLength 2000 \
    --numberOfProcessors 20 \
    --skipZeros -o ${data_path}${tissue}/${antibody}/matrix/${antibody}_higher_than_${bar_high}-W${window_size}-G${gap_size}-E${e_value}.mat.gz &

wait

plotProfile -m  ${data_path}${tissue}/${antibody}/matrix/${antibody}_less_than_${bar_low}-W${window_size}-G${gap_size}-E${e_value}.mat.gz \
            -out ${result_path}${tissue}/deeptools/${antibody}_less_than_${bar_low}-W${window_size}-G${gap_size}-E${e_value}-each-group.pdf --startLabel start --endLabel end \
            --numPlotsPerRow 5 --perGroup --legendLocation best \
            --plotTitle "${antibody}_less_than_${bar_low}" &
plotProfile -m  ${data_path}${tissue}/${antibody}/matrix/${antibody}_less_than_${bar_low}-W${window_size}-G${gap_size}-E${e_value}.mat.gz \
            -out ${result_path}${tissue}/deeptools/${antibody}_less_than_${bar_low}-W${window_size}-G${gap_size}-E${e_value}.pdf --startLabel start --endLabel end \
            --numPlotsPerRow 5  --legendLocation best \
            --plotTitle "${antibody}_less_than_${bar_low}" &

plotProfile -m  ${data_path}${tissue}/${antibody}/matrix/${antibody}_higher_than_${bar_high}-W${window_size}-G${gap_size}-E${e_value}.mat.gz \
            -out ${result_path}${tissue}/deeptools/${antibody}_higher_than_${bar_high}-W${window_size}-G${gap_size}-E${e_value}-each-group.pdf --startLabel start --endLabel end \
            --numPlotsPerRow 5 --perGroup --legendLocation best \
            --plotTitle "${antibody}_higher_than_${bar_high}" &
plotProfile -m  ${data_path}${tissue}/${antibody}/matrix/${antibody}_higher_than_${bar_high}-W${window_size}-G${gap_size}-E${e_value}.mat.gz \
            -out ${result_path}${tissue}/deeptools/${antibody}_higher_than_${bar_high}-W${window_size}-G${gap_size}-E${e_value}.pdf --startLabel start --endLabel end \
            --numPlotsPerRow 5  --legendLocation best \
            --plotTitle "${antibody}_higher_than_${bar_high}" &
