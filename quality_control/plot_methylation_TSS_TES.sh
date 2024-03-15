data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/20231201_LLX_Sample/bigWig/
result_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/20231201_LLX_Sample/
computeMatrix scale-regions -S ${data_path}LLX145*.bw ${data_path}LLX151*.bw\
  -R /storage/zhangyanxiaoLab/zhangyanxiao/annotations/hg38/gencode.v26.annotation.transcripts.tss1k.bed \
  --beforeRegionStartLength 2000 \
  --regionBodyLength 2000 \
  --afterRegionStartLength 2000 \
  --numberOfProcessors 20 \
  --skipZeros -o ${result_path}matrix/H3K9me3.mat.gz

plotProfile -m ${result_path}matrix/H3K9me3.mat.gz \
              -out ${result_path}heatmap/H3K9me3.pdf \
              --numPlotsPerRow 5 \
              --plotTitle "H3K9me3"

computeMatrix scale-regions -S ${data_path}LLX146*.bw ${data_path}LLX152*.bw\
  -R /storage/zhangyanxiaoLab/zhangyanxiao/annotations/hg38/gencode.v26.annotation.transcripts.tss1k.bed \
  --beforeRegionStartLength 2000 \
  --regionBodyLength 2000 \
  --afterRegionStartLength 2000 \
  --numberOfProcessors 20 \
  --skipZeros -o ${result_path}matrix/H3K36me3.mat.gz

plotProfile -m ${result_path}matrix/H3K36me3.mat.gz \
              -out ${result_path}heatmap/H3K36me3.pdf \
              --numPlotsPerRow 5 \
              --plotTitle "H3K36me3"

computeMatrix scale-regions -S ${data_path}LLX148*.bw ${data_path}LLX154*.bw\
  -R /storage/zhangyanxiaoLab/zhangyanxiao/annotations/hg38/gencode.v26.annotation.transcripts.tss1k.bed \
  --beforeRegionStartLength 2000 \
  --regionBodyLength 2000 \
  --afterRegionStartLength 2000 \
  --numberOfProcessors 20 \
  --skipZeros -o ${result_path}matrix/H3K27me3.mat.gz

plotProfile -m ${result_path}matrix/H3K27me3.mat.gz \
              -out ${result_path}heatmap/H3K27me3.pdf \
              --numPlotsPerRow 5 \
              --plotTitle "H3K27me3"
