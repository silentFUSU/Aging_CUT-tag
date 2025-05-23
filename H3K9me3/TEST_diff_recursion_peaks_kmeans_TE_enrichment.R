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
kmeans <- paste0("kmeans",c(1:4))
TE <- fread("~/ref_data/TE_reference/mm10_TE_metadata_20250116.bed")
setDT(TE)
setkey(TE,V1,V2,V3)
gtf <- import("~/ref_data/TE_reference/mm10_rmsk_TE.gtf", format = "gtf")
family_data <- as.data.frame(gtf[,c("gene_id","transcript_id","family_id","class_id")])
summary <- data.frame()
for(kmean in kmeans){
  df <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/",kmean,"_uinon_recursion_peaks.bed"))
  df$label <- paste0(df$V1,":",df$V2,"-",df$V3)
  length <- data.frame(label=df$label,length=df$V3-df$V2+1)
  df <- as.data.table(df)
  setDT(df)
  setkey(df,V1,V2,V3)
  overlaps <- foverlaps(df, TE, type = "any", nomatch = 0L)  
  overlaps <- as.data.frame(overlaps)
  overlaps$TE_length <- overlaps$V3 - overlaps$V2 +1
  result <- overlaps %>%
    group_by(label) %>%
    summarize(total_TE_length = sum(TE_length, na.rm = TRUE))
  
  t_summary <- as.data.frame(table(overlaps$label))
  t_summary$condition <- kmean
  colnames(t_summary)[1] <- "label"
  t_summary <- merge(t_summary,length,by="label")
  t_summary$enrichment <- t_summary$Freq/t_summary$length*1000
  t_summary <- merge(t_summary,result,by="label")
  t_summary$TE_length_percentage <- t_summary$total_TE_length/t_summary$length*100
  summary <- rbind(summary,t_summary[,c("condition","enrichment","TE_length_percentage")])  
}

ggplot(summary,aes(x=condition,y=enrichment,fill=condition))+    
  geom_boxplot()+
  ggtitle(paste0("H3K9me3 peaks TE enrichment"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("")+labs(fill = "", color = "") +ylab("Enrichment score")

ggplot(summary,aes(x=condition,y=TE_length_percentage,fill=condition))+    
  geom_boxplot()+
  ggtitle(paste0("H3K9me3 peaks TE length percentage"))+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("")+labs(fill = "", color = "") +ylab("TE length percentage")

classes <- c("LINE","LTR","SINE","DNA","RNA")
for(class in classes){
  summary <- data.frame()
  for(kmean in kmeans){
    df <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/",kmean,"_uinon_recursion_peaks.bed"))
    df <- df[which(df$V1 %in% paste0("chr",c(1:19,"X"))),]
    df$label <- paste0(df$V1,":",df$V2,"-",df$V3)
    length <- data.frame(label=df$label,length=df$V3-df$V2+1)
    df <- as.data.table(df)
    setDT(df)
    setkey(df,V1,V2,V3)
    overlaps <- foverlaps(df, TE, type = "any", nomatch = 0L)  
    overlaps <- overlaps[which(overlaps$V4 %in% family_data$transcript_id[which(family_data$class_id==class)]),]
    overlaps <- as.data.frame(overlaps)
    overlaps$TE_length <- overlaps$V3 - overlaps$V2 +1
    result <- overlaps %>%
      group_by(label) %>%
      summarize(total_TE_length = sum(TE_length, na.rm = TRUE))
    
    t_summary <- as.data.frame(table(overlaps$label))
    t_summary$condition <- kmean
    colnames(t_summary)[1] <- "label"
    t_summary <- merge(t_summary,length,by="label")
    t_summary$enrichment <- t_summary$Freq/t_summary$length*1000
    t_summary <- merge(t_summary,result,by="label")
    t_summary$TE_length_percentage <- t_summary$total_TE_length/t_summary$length*100
    summary <- rbind(summary,t_summary[,c("condition","enrichment","TE_length_percentage")])  
  }
  ggplot(summary,aes(x=condition,y=enrichment,fill=condition))+    
    geom_boxplot()+
    ggtitle(paste0("H3K9me3 peaks ",class," TE enrichment"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("")+labs(fill = "", color = "") +ylab("Enrichment score")
  
  ggplot(summary,aes(x=condition,y=TE_length_percentage,fill=condition))+    
    geom_boxplot()+
    ggtitle(paste0("H3K9me3 peaks ",class," TE length percentage"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("")+labs(fill = "", color = "") +ylab("TE length percentage")
}

families <- unique(family_data$family_id[which(family_data$class_id=="LTR")])
classes <- c("LINE","LTR","SINE","DNA","RNA")
for(family in families){
  summary <- data.frame()
  for(kmean in kmeans){
    df <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/",kmean,"_uinon_recursion_peaks.bed"))
    df <- df[which(df$V1 %in% paste0("chr",c(1:19,"X"))),]
    df$label <- paste0(df$V1,":",df$V2,"-",df$V3)
    length <- data.frame(label=df$label,length=df$V3-df$V2+1)
    df <- as.data.table(df)
    setDT(df)
    setkey(df,V1,V2,V3)
    overlaps <- foverlaps(df, TE, type = "any", nomatch = 0L)  
    overlaps <- overlaps[which(overlaps$V4 %in% family_data$transcript_id[which(family_data$family_id==family)]),]
    overlaps <- as.data.frame(overlaps)
    overlaps$TE_length <- overlaps$V3 - overlaps$V2 +1
    result <- overlaps %>%
      group_by(label) %>%
      summarize(total_TE_length = sum(TE_length, na.rm = TRUE))
    
    t_summary <- as.data.frame(table(overlaps$label))
    t_summary$condition <- kmean
    colnames(t_summary)[1] <- "label"
    t_summary <- merge(t_summary,length,by="label")
    t_summary$enrichment <- t_summary$Freq/t_summary$length*1000
    t_summary <- merge(t_summary,result,by="label")
    t_summary$TE_length_percentage <- t_summary$total_TE_length/t_summary$length*100
    summary <- rbind(summary,t_summary[,c("condition","enrichment","TE_length_percentage")])  
  }
  ggplot(summary,aes(x=condition,y=enrichment,fill=condition))+    
    geom_boxplot()+
    ggtitle(paste0("H3K9me3 peaks ",family," TE enrichment"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("")+labs(fill = "", color = "") +ylab("Enrichment score")
  
  ggplot(summary,aes(x=condition,y=TE_length_percentage,fill=condition))+    
    geom_boxplot()+
    ggtitle(paste0("H3K9me3 peaks ",family," TE length percentage"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("")+labs(fill = "", color = "") +ylab("TE length percentage")
}

for(family in c("ERVK","ERV1")){
  # df <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/kmeans1_uinon_recursion_peaks.bed"))
  # df <- df[which(df$V1 %in% paste0("chr",c(1:19,"X"))),]
  # df$label <- paste0(df$V1,":",df$V2,"-",df$V3)
  # length <- data.frame(label=df$label,length=df$V3-df$V2+1)
  # df <- as.data.table(df)
  # setDT(df)
  # setkey(df,V1,V2,V3)
  # overlaps <- foverlaps(df, TE, type = "any", nomatch = 0L)  
  # overlaps <- overlaps[which(overlaps$V4 %in% family_data$transcript_id[which(family_data$family_id==family)]),]
  # overlaps <- as.data.frame(overlaps)
  # overlaps <- merge(overlaps,family_data[,c("gene_id","transcript_id")],by.x="V4",by.y="transcript_id")
  # overlaps$TE_length <- overlaps$V3 - overlaps$V2 +1
  # result <- overlaps %>%
  #   group_by(label,gene_id) %>%
  #   summarize(total_TE_length = sum(TE_length, na.rm = TRUE))
  # t_summary <- result
  # t_summary <- merge(t_summary,length,by="label")
  # t_summary$TE_length_percentage <- t_summary$total_TE_length/t_summary$length*100
  # summary <- t_summary
  # gene_id_table <- as.data.frame(table(summary$gene_id))
  # summary <- summary[which(summary$gene_id %in% gene_id_table$Var1[which(gene_id_table$Freq>20)]),]
  # summary_median <- summary %>%
  #   group_by(gene_id) %>%
  #   summarise(median_TE_length_percentage = median(TE_length_percentage))
  # summary_median <- summary_median[order(summary_median$median_TE_length_percentage,decreasing = T),]
  # summary$gene_id <- factor(summary$gene_id,levels=summary_median$gene_id)
  # top_te_length_percentage <- summary_median$gene_id[1:10]
  # ggplot(summary[which(summary$gene_id %in% summary_median$gene_id[1:10]),],aes(x=gene_id,y=TE_length_percentage,fill=gene_id))+    
  #   geom_boxplot()+
  #   ggtitle(paste0("H3K9me3 peaks ",family," TE length percentage"))+
  #   theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")+
  #   xlab("")+labs(fill = "", color = "") +ylab("TE length percentage")
  # 
  # result <- overlaps %>%
  #   group_by(label, gene_id) %>%
  #   summarise(count = n(), .groups = "drop")
  # t_summary <- result
  # t_summary <- merge(t_summary,length,by="label")
  # t_summary$enrichment <- t_summary$count/t_summary$length*1000
  # summary <- t_summary
  # summary_median <- summary %>%
  #   group_by(gene_id) %>%
  #   summarise(enrichment = median(enrichment))
  # summary_median <- summary_median[order(summary_median$enrichment,decreasing = T),]
  # summary$gene_id <- factor(summary$gene_id,levels=summary_median$gene_id)
  # top_enrichment <- summary_median$gene_id[1:10]
  # ggplot(summary[which(summary$gene_id %in% summary_median$gene_id[1:10]),],aes(x=gene_id,y=enrichment,fill=gene_id))+    
  #   geom_boxplot()+
  #   ggtitle(paste0("H3K9me3 peaks ",family," TE enrichment"))+
  #   theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")+
  #   xlab("")+labs(fill = "", color = "") +ylab("Enrichment score")
  # 
  kmeans_summary <- data.frame()
  for(kmean in kmeans){
    df <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/",kmean,"_uinon_recursion_peaks.bed"))
    # df <- df[which(df$V1 %in% paste0("chr",c(1:19,"X"))),]
    df$label <- paste0(df$V1,":",df$V2,"-",df$V3)
    length <- data.frame(label=df$label,length=df$V3-df$V2+1)
    df <- as.data.table(df)
    setDT(df)
    setkey(df,V1,V2,V3)
    overlaps <- foverlaps(df, TE, type = "any", nomatch = 0L)  
    overlaps <- overlaps[which(overlaps$V4 %in% family_data$transcript_id[which(family_data$family_id==family)]),]
    overlaps <- as.data.frame(overlaps)
    overlaps <- merge(overlaps,family_data[,c("gene_id","transcript_id")],by.x="V4",by.y="transcript_id")
    overlaps$TE_length <- overlaps$V3 - overlaps$V2 +1
    result <- overlaps %>%
      group_by(label,gene_id) %>%
      summarize(total_TE_length = sum(TE_length, na.rm = TRUE))

    t_summary <- result
    t_summary$condition <- kmean
    colnames(t_summary)[1] <- "label"
    t_summary <- merge(t_summary,length,by="label")
    t_summary$TE_length_percentage <- t_summary$total_TE_length/t_summary$length*100
    kmeans_summary <- rbind(kmeans_summary,t_summary[,c("condition","gene_id","TE_length_percentage")])  
  }
  gene_id_table <- as.data.frame(table(kmeans_summary$gene_id))
  kmeans_summary <- kmeans_summary[which(kmeans_summary$gene_id %in% gene_id_table$Var1[which(gene_id_table$Freq>20)]),]
  kmeans_summary_median <- kmeans_summary %>%
    group_by(condition,gene_id) %>%
    summarise(median_TE_length_percentage = median(TE_length_percentage))
  to_plot <- reshape2::dcast(kmeans_summary_median, gene_id ~ condition, value.var = "median_TE_length_percentage")
  rownames(to_plot) <- to_plot$gene_id
  to_plot <- to_plot[,-1]
  breaks <- c(seq(0,0.5, length.out = 40))
  pheatmap::pheatmap(to_plot,main = paste0(family," family"),breaks = breaks)
  # kmeans_summary$gene_id <- factor(kmeans_summary$gene_id,levels=top_te_length_percentage)
  # ggplot(kmeans_summary, aes(x=condition,y=TE_length_percentage,fill=gene_id))+    
  #   geom_boxplot()+
  #   ggtitle(paste0("H3K9me3 peaks ",family," TE length percentage"))+
  #   theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  #   xlab("")+labs(fill = "", color = "") +ylab("TE length percentage")
  

  kmeans_summary <- data.frame()
  for(kmean in kmeans){
    df <- read.table(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/",kmean,"_uinon_recursion_peaks.bed"))
    # df <- df[which(df$V1 %in% paste0("chr",c(1:19,"X"))),]
    df$label <- paste0(df$V1,":",df$V2,"-",df$V3)
    length <- data.frame(label=df$label,length=df$V3-df$V2+1)
    df <- as.data.table(df)
    setDT(df)
    setkey(df,V1,V2,V3)
    overlaps <- foverlaps(df, TE, type = "any", nomatch = 0L)  
    overlaps <- overlaps[which(overlaps$V4 %in% family_data$transcript_id[which(family_data$family_id==family )]),]
    overlaps <- as.data.frame(overlaps)
    overlaps <- merge(overlaps,family_data[,c("gene_id","transcript_id")],by.x="V4",by.y="transcript_id")
    overlaps$TE_length <- overlaps$V3 - overlaps$V2 +1
    result <- overlaps %>%
      group_by(label, gene_id) %>%
      summarise(count = n(), .groups = "drop")
    
    t_summary <- result
    t_summary$condition <- kmean
    colnames(t_summary)[1] <- "label"
    t_summary <- merge(t_summary,length,by="label")
    t_summary$enrichment <- t_summary$count/t_summary$length*1000
    kmeans_summary <- rbind(kmeans_summary,t_summary[,c("condition","gene_id","enrichment")])  
  }
  gene_id_table <- as.data.frame(table(kmeans_summary$gene_id))
  kmeans_summary <- kmeans_summary[which(kmeans_summary$gene_id %in% gene_id_table$Var1[which(gene_id_table$Freq>20)]),]
  kmeans_summary_median <- kmeans_summary %>%
    group_by(condition,gene_id) %>%
    summarise(median_enrichment = median(enrichment))
  to_plot <- reshape2::dcast(kmeans_summary_median, gene_id ~ condition, value.var = "median_enrichment")
  rownames(to_plot) <- to_plot$gene_id
  to_plot <- to_plot[,-1]
  # breaks <- c(seq(0,0.5, length.out = 40))
  pheatmap::pheatmap(to_plot,main = paste0(family," family"),scale="row")
  # kmeans_summary$gene_id <- factor(kmeans_summary$gene_id,levels=top_enrichment)
  # ggplot(kmeans_summary, aes(x=condition,y=enrichment,fill=gene_id))+    
  #   geom_boxplot()+
  #   ggtitle(paste0("H3K9me3 peaks ",family," TE enrichment"))+
  #   theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  #   xlab("")+labs(fill = "", color = "") +ylab("Enrichment score")
}

