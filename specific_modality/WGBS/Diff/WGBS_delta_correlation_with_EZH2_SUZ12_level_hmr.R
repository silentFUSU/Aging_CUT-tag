rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
library(ggsignif)
library(data.table)
library(dplyr)
library(GenomeInfoDb)
library("GenomicRanges")
library(genomation)
tissue_label_change <- function(tissue){
  if(tissue=="brain"){
    tissue_label <- "Cortex"
  }else if(tissue == "Hip"){
    tissue_label <- "Hippocampus"
  }else if(tissue == "CB"){
    tissue_label <- "Cerebellum"
  }else{
    tissue_label <- str_to_title(tissue)
    if(tissue_label == "Bonemarrow"){
      tissue_label <- "Bone Marrow"
    }else if(tissue_label == "Bat"){
      tissue_label <- "BAT"
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}

tissue <- "stomach"
hmr <- read.table(paste0("data/samples/WGBS/",tissue,"/hmr/all_samples_hmr.bed"))
if(tissue %in% c("mammarygland","ovary","uterus")){
  hmr <- hmr[which(hmr$V1 %in% paste0("chr",c(1:19,"X"))),]
}else{
  hmr <- hmr[which(hmr$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
}
blacklist <- read.table("~/ref_data/mm10-blacklist.v2.bed",sep = "\t")
blacklist <- as.data.table(blacklist)
setDT(blacklist)
setkey(blacklist,V1,V2,V3)

hmr_regions <- as.data.table(hmr)
setDT(hmr_regions)
setkey(hmr_regions,V1,V2,V3)
overlaps <- foverlaps(hmr_regions, blacklist, type = "any", nomatch = 0L)  

hmr <- hmr[which(!hmr$V4 %in% overlaps$i.V4),]

EZH2 <- fread("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_input/EZH2_comp.bedGraph")
colnames(EZH2)[4] <- "EZH2"

SUZ12 <- fread("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_input/SUZ12_comp.bedGraph")
colnames(SUZ12)[4] <- "SUZ12"


hmr_regions <- as.data.table(hmr)
setDT(hmr_regions)
setkey(hmr_regions,V1,V2,V3)

setDT(EZH2)
setkey(EZH2,V1,V2,V3)


setDT(SUZ12)
setkey(SUZ12,V1,V2,V3)

overlaps <- foverlaps(EZH2,hmr_regions, type = "any", nomatch = 0L) 
overlaps<- overlaps[, .(EZH2 = mean(EZH2)), by = .(V1, V2, V3, V4)]
overlaps <- as.data.table(overlaps[,c("V1","V2","V3","V4","EZH2")])

setDT(overlaps)
setkey(overlaps,V1,V2,V3)

overlaps_SUZ12 <- foverlaps(SUZ12,overlaps, type = "any", nomatch = 0L)  
overlaps_SUZ12<- overlaps_SUZ12[, .(SUZ12 = mean(SUZ12)), by = .(V1, V2, V3, V4, EZH2)]

overlaps_SUZ12$score <- (overlaps_SUZ12$SUZ12 + overlaps_SUZ12$EZH2)/2

overlaps_SUZ12 <- as.data.frame(overlaps_SUZ12)
overlaps_SUZ12 <- overlaps_SUZ12[,c("V1","V2","V3","V4","score")]

score <- overlaps_SUZ12[,c("V4","score")]

regions <- as.data.table(hmr[,c(1:4)])
setDT(regions)
setkey(regions,V1,V2,V3)
search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
search_table <- search_table[which(search_table$tissue==tissue),]
summary <- data.frame()
for(sample in search_table$sample_name){
  df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
  setDT(df)
  setkey(df,V1,V2,V3)  
  overlaps <- foverlaps(df, regions, type = "any", nomatch = 0L)  
  result <- overlaps[, .(V4_sum = sum(i.V4), V5_sum = sum(V5)),by=V4]
  result <- as.data.frame(result)
  result$methylation <- result$V4_sum/result$V5_sum*100
  result <- result[which(result$V5_sum >=20),]
  result <- result[,c("V4","methylation")]
  colnames(result) <- c("LMR",sample)
  if(nrow(summary)==0){
    summary <- result
  }else{
    summary <- merge(summary,result,by="LMR")
  }
}

young <- summary[,c("LMR",search_table$sample_name[which(search_table$age=="3M")])]
young$young <- rowMeans(young[,-1])

old <- summary[,c("LMR",search_table$sample_name[which(search_table$age=="24M")])]
old$old <- rowMeans(old[,-1])

methylation <- merge(young[,c("LMR","young")],old[,c("LMR","old")],by="LMR")
methylation$delta <- methylation$old - methylation$young

to_plot <- merge(methylation,score,by.x="LMR",by.y="V4")
to_plot$condition <- "all"
to_plot <- to_plot[order(to_plot$score,decreasing = T),]
to_plot$condition[1:1000] <- "top1000"
color <- setNames(c("blue","red"),c("all","top1000"))
p <- ggplot(to_plot, aes(x = score, y = delta, color=condition)) +
  geom_point(size=0.5) +
  scale_color_manual(values=color)+
  labs(title = tissue_label_change(tissue),
       x = "PRC2 binding level (Fold change)",
       y = "Delta") +
  theme_bw() +
  scale_x_log10(
    limits = c(0.9, 120), 
    breaks = c(1, 10, 100),
    labels = c("10^0", "10^1", "10^2")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "green")
ggsave("result/figures/Stomach_WGBS_delta_correlation_with_EZH2_SUZ12.pdf",p,width = 7,height = 6)
### find top1000
hmr <- read.table(paste0("data/samples/WGBS/",tissue,"/hmr/all_samples_hmr.bed"))
if(tissue %in% c("mammarygland","ovary","uterus")){
  hmr <- hmr[which(hmr$V1 %in% paste0("chr",c(1:19,"X"))),]
}else{
  hmr <- hmr[which(hmr$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
}
blacklist <- read.table("~/ref_data/mm10-blacklist.v2.bed",sep = "\t")
blacklist <- as.data.table(blacklist)
setDT(blacklist)
setkey(blacklist,V1,V2,V3)

hmr_regions <- as.data.table(hmr)
setDT(hmr_regions)

setkey(hmr_regions,V1,V2,V3)
overlaps <- foverlaps(hmr_regions, blacklist, type = "any", nomatch = 0L)  

hmr <- hmr[which(!hmr$V4 %in% overlaps$i.V4),]

EZH2 <- read.table("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_input/EZH2_peaks.narrowPeak")

EZH2 <- EZH2[,c(1:3,7,9)]
colnames(EZH2)[4] <- "EZH2_signal"
colnames(EZH2)[5] <- "EZH2"

SUZ12 <- read.table("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_input/SUZ12_peaks.narrowPeak")
SUZ12 <- SUZ12[,c(1:3,7,9)]
colnames(SUZ12)[4] <- "SUZ12_signal"
colnames(SUZ12)[5] <- "SUZ12"

hmr_regions <- as.data.table(hmr)
setDT(hmr_regions)
setkey(hmr_regions,V1,V2,V3)

EZH2 <- as.data.table(EZH2)
setDT(EZH2)
setkey(EZH2,V1,V2,V3)

SUZ12 <- as.data.table(SUZ12)
setDT(SUZ12)
setkey(SUZ12,V1,V2,V3)

overlaps <- foverlaps(EZH2,hmr_regions, type = "any", nomatch = 0L)  
overlaps <- as.data.table(overlaps[,c("V1","V2","V3","V4","EZH2_signal","EZH2")])
setDT(overlaps)
setkey(overlaps,V1,V2,V3)

overlaps_SUZ12 <- foverlaps(overlaps, SUZ12, type = "any", nomatch = 0L)  
overlaps_SUZ12$SUZ12 <- as.numeric(overlaps_SUZ12$SUZ12)
overlaps_SUZ12$EZH2 <- as.numeric(overlaps_SUZ12$EZH2)

overlaps_SUZ12$score <- (overlaps_SUZ12$SUZ12 + overlaps_SUZ12$EZH2)/2
overlaps_SUZ12$enrichment <- (overlaps_SUZ12$SUZ12_signal + overlaps_SUZ12$EZH2_signal)/2

overlaps_SUZ12 <- as.data.frame(overlaps_SUZ12)

score <- overlaps_SUZ12 %>%
  group_by(V4) %>%
  summarize(
    score_max = max(score, na.rm = TRUE),
    enrichment_max = max(enrichment, na.rm = TRUE)
  )
top1000 <- score[order(score$score_max,decreasing = T),]
top1000 <- top1000$V4[1:1000]


to_plot$condition <- "all"
to_plot$condition[which(to_plot$LMR %in% top1000)] <- "top1000"
color <- setNames(c("blue","red"),c("all","top1000"))
ggplot(to_plot, aes(x = score, y = delta, color=condition)) +
  geom_point(size=0.5) +
  scale_color_manual(values=color)+
  labs(title = tissue_label_change(tissue),
       x = "PRC2 binding level (Fold change)",
       y = "Delta") +
  theme_bw() +
  scale_x_log10(
    limits = c(0.9, 120), 
    breaks = c(1, 10, 100),
    labels = c("10^0", "10^1", "10^2")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "green")




