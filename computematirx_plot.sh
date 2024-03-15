data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/
tissue=brain
base_antibody=H3K27me3
antibody=H3K9me3
ref=mm10
window_size=1000
gap_size=3000
e_value=100
condition=up
bw=$(ls ${data_path}${tissue}/${antibody}/bw/*nodup.bw)

computeMatrix scale-regions -S ${bw} \
  -R ${data_path}${tissue}/${base_antibody}/bed/${base_antibody}_merge-W${window_size}-G${gap_size}-E${e_value}_${condition}.bed     \
  --beforeRegionStartLength 2000 --startLabel start --endLabel end \
  --regionBodyLength 2000 \
  --afterRegionStartLength 2000 \
  --numberOfProcessors 10 \
  --skipZeros -o ${data_path}${tissue}/${antibody}/matrix/${antibody}_base${base_antibody}_W${window_size}-G${gap_size}-E${e_value}_${condition}.mat.gz 

plotProfile -m  ${data_path}${tissue}/${antibody}/matrix/${antibody}_base${base_antibody}_W${window_size}-G${gap_size}-E${e_value}_${condition}.mat.gz    \
              -out ${result_path}${tissue}/deeptools/${antibody}_base${base_antibody}_W${window_size}-G${gap_size}-E${e_value}_${condition}.pdf --startLabel start --endLabel end \
              --numPlotsPerRow 5 --perGroup --legendLocation best \
              --plotTitle "${antibody}"

bw=$(ls ${data_path}${tissue}/${base_antibody}/bw/*nodup.bw)

computeMatrix scale-regions -S ${bw} \
  -R ${data_path}${tissue}/${base_antibody}/bed/${base_antibody}_merge-W${window_size}-G${gap_size}-E${e_value}_${condition}.bed     \
  --beforeRegionStartLength 2000 --startLabel start --endLabel end \
  --regionBodyLength 2000 \
  --afterRegionStartLength 2000 \
  --numberOfProcessors 10 \
  --skipZeros -o ${data_path}${tissue}/${base_antibody}/matrix/${base_antibody}_base${base_antibody}_W${window_size}-G${gap_size}-E${e_value}_${condition}.mat.gz 

plotProfile -m  ${data_path}${tissue}/${base_antibody}/matrix/${base_antibody}_base${base_antibody}_W${window_size}-G${gap_size}-E${e_value}_${condition}.mat.gz    \
              -out ${result_path}${tissue}/deeptools/${base_antibody}_base${base_antibody}_W${window_size}-G${gap_size}-E${e_value}_${condition}.pdf --startLabel start --endLabel end \
              --numPlotsPerRow 5 --perGroup --legendLocation best \
              --plotTitle "${base_antibody}"
