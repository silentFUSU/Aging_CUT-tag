rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(limma)
library(data.table)
resolution <- "20000"
TAD <- read.csv("data/public_data/WANG_cellular_aging_HiC/TAD/insulation_score/redundant_TAD.csv")
TAD <- TAD[,c(1:3)]
TAD_regions <- data.frame(chr=as.character(),start=as.numeric(),end=as.numeric())
i=2
while(i <= nrow(TAD)){
  if(TAD[i,"chr"] == TAD[(i-1),"chr"]){
    t_TAD_regions <- data.frame(chr=TAD$chr[i],start=TAD$start[i-1],end=TAD$start[i]) 
    TAD_regions <- rbind(TAD_regions,t_TAD_regions)
    i <- i+1
  }else{
    i <- i+1
  }
}
TAD_regions$start <- TAD_regions$start+1
TAD_regions <- as.data.table(TAD_regions)
setDT(TAD_regions)
setkey(TAD_regions,chr,start,end)

bed <- read.table(paste0("data/public_data/WANG_cellular_aging_HiC/raw_matrix/G1_",resolution,"_abs.bed"))
bed <- bed[which(bed$V1 %in% paste0("chr",c(1:22,"X"))),]
bed$V5 <- "NA"
count <- 1
bed[1,"V5"] <- count
for(i in c(2:nrow(bed))){
  if(bed[i,"V1"] != bed[i-1,"V1"]){
    count <- 1
    bed[i,"V5"] <- count
  }else{
    count <- count+1
    bed[i,"V5"] <- count
  }
}
bed$V2 <- bed$V2+1
bed <- as.data.table(bed)
setDT(bed)
setkey(bed,V1,V2,V3)

overlaps <- foverlaps(TAD_regions, bed, type = "any", nomatch = 0L)  
overlaps$TAD_label <- paste(overlaps$chr,overlaps$start,overlaps$end,sep = "-")
re <-  read.table(paste0("data/public_data/WANG_cellular_aging_HiC/differential_analysis/G_DS_20000.FDR"))
re$Significant <- "Stable"
re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 < 0)] <- "Down"
re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 > 0)] <- "Up"
re_sig <- re[which(re$Significant!="Stable"),]

# re_sig <- re_sig[which(abs(re_sig$V2 - re_sig$V3)>4),]
re_sig$V1 <- paste0("chr",re_sig$V1)
re_sig$V1[which(re_sig$V1=="chr23")] <- "chrX"
# re_sig_reverse <- data.frame(V1=re_sig$V1,V2=re_sig$V3, V3=re_sig$V2, V4=re_sig$V4,V5=re_sig$V5,V6=re_sig$V6,Significant=re_sig$Significant)
# re_sig <- rbind(re_sig,re_sig_reverse)
re_sig$label1 <- paste(re_sig$V1,re_sig$V2,sep = "-")
re_sig$label2 <- paste(re_sig$V1,re_sig$V3,sep = "-")
overlaps$label <- paste(overlaps$chr,overlaps$V5,sep = "-")

re_sig <- merge(re_sig,overlaps[,c("label","TAD_label")],by.x="label1",by.y="label")
colnames(re_sig)[which(colnames(re_sig)=="TAD_label")] <- "TAD_label1"

re_sig <- merge(re_sig,overlaps[,c("label","TAD_label")],by.x="label2",by.y="label")
colnames(re_sig)[which(colnames(re_sig)=="TAD_label")] <- "TAD_label2"
re_sig <- re_sig[which(re_sig$TAD_label1==re_sig$TAD_label2),]
result <- re_sig %>%  
  group_by(TAD_label1) %>%  
  summarise(  
    Up = sum(Significant == "Up"),  
    Down = sum(Significant == "Down")  
  )  
result$condition <- "Stable"
result$condition[which(result$Up - result$Down >=5)] <- "Up"
result$condition[which(result$Down - result$Up >=5)] <- "Down"
color <- setNames(c("red","blue","grey"),c("Up","Down","Stable"))
ggplot(result, aes(x = Down, y = Up, color=condition)) +      
  geom_point(size = 1) +  
  scale_color_manual(values = color) +
  labs(title = "DS vs. G",  
       x = "# of increased interactions in TADs",  
       y = "# of decreased interactions in TADs") +        
  xlim(0,80) +
  ylim(0,125) +
  theme_minimal()+
  theme(  
    # 修改所有字体的大小  
    plot.title = element_text(size = 16),  # 图标题字体大小  
    axis.title = element_text(size = 14),  # 轴标题字体大小  
    axis.text = element_text(size = 12),   # 坐标刻度标签字体大小  
    legend.text = element_text(size = 12), # 图例字体大小  
    legend.title = element_text(size = 14) # 图例标题字体大小  
  )  
  
antibody <- "H3K27me3"
histone <- read.table(paste0("data/public_data/cellular_aging_GSE133292/Chip_seq/",antibody,"_10kb_bins.counts"),header = T)
counts <- histone[,c(7:ncol(histone))]
rownames(counts)= histone$Geneid
pattern <- ".*bam\\.(SRR[0-9]+).*"
colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
lengths_kb <- histone$Length/1000
rpkm_matrix <- matrix(0, nrow = nrow(counts), ncol = ncol(counts))  
rownames(rpkm_matrix) <- rownames(counts)  
colnames(rpkm_matrix) <- colnames(counts) 
total_mapped_reads <- colSums(counts)
for (i in 1:ncol(counts)) {  
  rpkm_matrix[, i] <- (counts[, i] / lengths_kb) / (total_mapped_reads[i] / 1e6)  
} 
rpkm_matrix <- as.data.frame(rpkm_matrix)
G_rpkm_matrix <- rpkm_matrix[,c(1:2)]
DS_rpkm_matrix <- rpkm_matrix[,c(3:4)]
G_rpkm_matrix <- as.data.frame(rowMeans(G_rpkm_matrix))
DS_rpkm_matrix <- as.data.frame(rowMeans(DS_rpkm_matrix))
rpkm_matrix <- cbind(G_rpkm_matrix,DS_rpkm_matrix)
colnames(rpkm_matrix) <- c("G","DS")
histone <- cbind(histone[,c(1:6)],rpkm_matrix)
histone$Start <- histone$Start+1
histone <- as.data.table(histone)
setDT(histone)
setkey(histone,Chr,Start,End)
overlaps_histone <- foverlaps(TAD_regions, histone, type = "any", nomatch = 0L) 
overlaps_histone$TAD_label <- paste(overlaps_histone$chr,overlaps_histone$start,overlaps_histone$end,sep = "-")

overlaps_histone <- overlaps_histone %>%  
  group_by(TAD_label) %>%  
  mutate(G_enrichment = median(G),
         DS_enrichment = median(DS)) %>%  
  ungroup() %>%  
  select(TAD_label, G_enrichment, DS_enrichment) %>%  
  distinct(TAD_label, .keep_all = TRUE) 

to_plot <- merge(result,overlaps_histone,by.x="TAD_label1", by.y="TAD_label")
ggplot(to_plot, aes(x = G_enrichment, y = DS_enrichment, color=condition)) +      
  geom_point(size = 0.5) +  
  scale_color_manual(values = color) +
  labs(title = antibody,  
       x = "Median RPKM in G",  
       y = "Median RPKM in DS",
       color = "TAD types") +        
  theme_minimal()+
  xlim(0,1.5) +
  ylim(0,1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "green") +
  theme(  
    # 修改所有字体的大小  
    plot.title = element_text(size = 16),  # 图标题字体大小  
    axis.title = element_text(size = 14),  # 轴标题字体大小  
    axis.text = element_text(size = 12),   # 坐标刻度标签字体大小  
    legend.text = element_text(size = 12), # 图例字体大小  
    legend.title = element_text(size = 14) # 图例标题字体大小  
  )  
