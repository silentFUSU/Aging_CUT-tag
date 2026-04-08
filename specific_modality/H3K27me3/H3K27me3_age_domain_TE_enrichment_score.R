rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(rtracklayer)
options(bitmapType="cairo")  

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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}

TE <- fread("~/ref_data/TE_reference/mm10_TE_metadata_20250116.bed")
setDT(TE)
setkey(TE,V1,V2,V3)
gtf <- import("~/ref_data/TE_reference/mm10_rmsk_TE.gtf", format = "gtf")
family_data <- as.data.frame(gtf[,c("gene_id","transcript_id","family_id","class_id")])
tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
Classification <- "family_id"
summary_enrichment <- data.frame()

for(tissue in tissues){
  df <- read.table(paste0("data/samples/",tissue,"/H3K27me3/peaks/edd/age_domain.bed"))
  colnames(df)[4] <- "label"
  df$V2 <- df$V2 + 1
  length <- data.frame(label=df$label,length=df$V3-df$V2+1)
  df <- as.data.table(df)
  setDT(df)
  setkey(df,V1,V2,V3)
  overlaps <- foverlaps(df, TE, type = "any", nomatch = 0L)  
  overlaps <- merge(overlaps,family_data[,c("transcript_id","gene_id","family_id","class_id")],by.x="V4",by.y="transcript_id")
  
  overlaps <- as.data.frame(overlaps)
  overlaps$TE_length <- overlaps$V3 - overlaps$V2 +1
  colnames(overlaps)[which(colnames(overlaps)==Classification)] <- "classification"
  overlaps$classification <- paste0(overlaps$classification,"-",overlaps$class_id)
  
  t_summary_enrichment <- overlaps %>%
    group_by(label, classification) %>%
    summarise(count = n(), .groups = "drop")
  
  t_summary_enrichment <- reshape2::dcast(t_summary_enrichment, classification ~ label, value.var = "count")
  t_summary_enrichment[is.na(t_summary_enrichment)] <- 0
  t_summary_enrichment <- reshape2::melt(t_summary_enrichment)
  
  colnames(t_summary_enrichment)[2] <- "label"
  colnames(t_summary_enrichment)[3] <- "count"
  t_summary_enrichment$tissue <- tissue_label_change(tissue) 
  t_summary_enrichment <- merge(t_summary_enrichment,length,by="label")
  t_summary_enrichment$enrichment <- t_summary_enrichment$count/t_summary_enrichment$length*1000
  
  t_summary_enrichment_mean <- t_summary_enrichment %>%
    group_by(classification) %>%
    summarise(mean_enrichment = mean(enrichment))
  
  colnames(t_summary_enrichment_mean)[2] <- tissue_label_change(tissue)
  
  if(nrow(summary_enrichment) == 0){
    summary_enrichment <- t_summary_enrichment_mean
  }else{
    summary_enrichment <- merge(summary_enrichment,t_summary_enrichment_mean,by="classification")
  }
}
to_plot <- summary_enrichment
rownames(to_plot) <- summary_enrichment$classification
to_plot <- to_plot[,-1]

pheatmap::pheatmap(to_plot,scale = "column")

#### L1
family_data <- family_data[which(family_data$family_id %in% c("L1")),]
Classification <- "gene_id"
summary_enrichment <- data.frame()
for(tissue in tissues){
  df <- read.table(paste0("data/samples/",tissue,"/H3K27me3/peaks/edd/age_domain.bed"))
  df$label <- paste0(df$V1,":",df$V2,"-",df$V3)
  length <- data.frame(label=df$label,length=df$V3-df$V2+1)
  df <- as.data.table(df)
  setDT(df)
  setkey(df,V1,V2,V3)
  
  overlaps <- foverlaps(df, TE, type = "any", nomatch = 0L)  
  overlaps <- merge(overlaps,family_data[,c("transcript_id","gene_id","family_id","class_id")],by.x="V4",by.y="transcript_id")
  
  overlaps <- as.data.frame(overlaps)
  overlaps$TE_length <- overlaps$V3 - overlaps$V2 +1
  colnames(overlaps)[which(colnames(overlaps)==Classification)] <- "classification"
  overlaps$classification <- paste0(overlaps$classification,"-",overlaps$class_id)
  
  t_summary_enrichment <- overlaps %>%
    group_by(label, classification) %>%
    summarise(count = n(), .groups = "drop")
  
  t_summary_enrichment <- reshape2::dcast(t_summary_enrichment, classification ~ label, value.var = "count")
  t_summary_enrichment[is.na(t_summary_enrichment)] <- 0
  t_summary_enrichment <- reshape2::melt(t_summary_enrichment)
  
  colnames(t_summary_enrichment)[2] <- "label"
  colnames(t_summary_enrichment)[3] <- "count"
  t_summary_enrichment$tissue <- tissue_label_change(tissue) 
  t_summary_enrichment <- merge(t_summary_enrichment,length,by="label")
  t_summary_enrichment$enrichment <- t_summary_enrichment$count/t_summary_enrichment$length*1000
  
  t_summary_enrichment_mean <- t_summary_enrichment %>%
    group_by(classification) %>%
    summarise(mean_enrichment = mean(enrichment))
  
  colnames(t_summary_enrichment_mean)[2] <- tissue_label_change(tissue)
  
  if(nrow(summary_enrichment) == 0){
    summary_enrichment <- t_summary_enrichment_mean
  }else{
    summary_enrichment <- merge(summary_enrichment,t_summary_enrichment_mean,by="classification")
  }
}
to_plot <- summary_enrichment
rownames(to_plot) <- summary_enrichment$classification
to_plot <- to_plot[,-1]
pheatmap::pheatmap(to_plot)
