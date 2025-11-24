rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
library(clusterProfiler)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(tidyverse)
library(data.table)
library(GO.db)
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
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissue_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv")) 

  bin <- as.data.table(df[,c(1:4)])
  setDT(bin)
  setkey(bin,Chr,Start,End)
  peaks <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
  peaks <- as.data.table(peaks)
  setDT(peaks)
  setkey(peaks,V1,V2,V3)
  overlaps <- foverlaps(bin, peaks, type = "any", nomatch = 0L)  
  df <- df[which(df$Significant == "Down" & df$Geneid %in% overlaps$Geneid),]
  bin <- as.data.table(df[,c(1:4)])
  setDT(bin)
  setkey(bin,Chr,Start,End)
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue),]
  CG_summary <- data.frame()
  for(sample in search_table$sample_name){
    CG <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    setDT(CG)
    setkey(CG,V1,V2,V3)  
    overlaps <- foverlaps(CG, bin, type = "any", nomatch = 0L)      
    result <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(V5)), by = Geneid]
    result$methylation <- result$V4_sum/result$V5_sum *100
    result$sample <- sample
    result$age <- search_table$age[which(search_table$sample_name==sample)]
    CG_summary <- rbind(CG_summary,result)
  }
  CG_summary <- CG_summary %>%
    group_by(age, Geneid) %>%
    summarize(methylation = mean(methylation, na.rm = TRUE))
  CG_summary$tissue <- tissue_label_change(tissue)
  RNA_df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  CG_summary <- merge(CG_summary,RNA_df[,c("X","Significant")],by.x="Geneid",by.y="X")
  tissue_summary <- rbind(tissue_summary,CG_summary)
}
to_plot <- tissue_summary %>%
  group_by(Geneid, age, tissue, Significant) %>%
  summarize(avg_methylation = mean(methylation, na.rm = TRUE))

to_plot <- to_plot %>%
  group_by(age, tissue, Significant) %>%
  summarize(avg_methylation = mean(avg_methylation, na.rm = TRUE))
to_plot$age <- factor(to_plot$age, levels = c("3M","24M"))
to_plot$Significant <- factor(to_plot$Significant,levels = c("Up","Stable","Down"))
ggplot(to_plot, aes(x = Significant, y = avg_methylation,fill=age)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw()+  
  ylab("RPKM")+
  xlab(NULL)+
  scale_fill_brewer(palette = "Pastel1") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
    axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
    axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
    axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
    legend.text = element_text(size = 12)
  ) +ylim(0,100)
p_value_summary <- data.frame()
for(condition in c("Up","Stable","Down")){
  df <- to_plot[which(to_plot$Significant==condition),]
  df$label <- paste0(df$tissue,"-",df$Significant)
  # df$label <- paste0(df$tissue,"-",df$Geneid)
  df_young <- df[which(df$age=="3M"),]
  df_old <- df[which(df$age=="24M"),]
  df <- merge(df_young[,c("label","avg_methylation")],df_old[,c("label","avg_methylation")],by="label")
  test <- t.test(df$avg_methylation.x,df$avg_methylation.y,paired = T)
  t_p_value_summary <- data.frame(condition=condition,p_value=test$p.value)
  p_value_summary <- rbind(p_value_summary,t_p_value_summary)
}
head(p_value_summary)

#overlap with CgG island
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
tissue_summary <- data.frame()
cpg <- read.table("~/ref_data/for_normal_mapping/mm10/cpgi.mm10.bed.txt")
cpg <- cpg[which(cpg$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
cpg <- as.data.table(cpg[,c(1:3)])
setDT(cpg)
setkey(cpg,V1,V2,V3)
genes <- read.table("~/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS.bed")
genes <- as.data.table(genes[,c(1:3,6)])
setDT(genes)
setkey(genes,V1,V2,V3)
CpG_overlaps <- foverlaps(genes, cpg, type = "any", nomatch = 0L)  
tissue_summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv")) 
  bin <- as.data.table(df[,c(1:4)])
  setDT(bin)
  setkey(bin,Chr,Start,End)
  peaks <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
  peaks <- as.data.table(peaks)
  setDT(peaks)
  setkey(peaks,V1,V2,V3)
  overlaps <- foverlaps(bin, peaks, type = "any", nomatch = 0L)  
  df <- df[which(df$Significant == "Down" & df$Geneid %in% overlaps$Geneid),]
  df$condition <- "out"
  df$condition[which(df$Geneid %in% CpG_overlaps$V6)] <- "overlap"
  df <- df[,c("Geneid","condition")]
  RNA_df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  df <- merge(df,RNA_df[,c("X","Significant")],by.x="Geneid",by.y="X")
  result <- df %>%
    group_by(Significant) %>%
    summarize(
      overlap_count = sum(condition == "overlap"),
      out_count = sum(condition == "out")
    )
  result$overlap <- result$overlap_count/(result$overlap_count+result$out_count)*100
  result$tissue <- tissue_label_change(tissue)
  tissue_summary <- rbind(tissue_summary,result)
}
tissue_summary$out <- 100 - tissue_summary$overlap
condition <- "Stable"
to_plot <- tissue_summary[which(tissue_summary$Significant==condition),]

to_plot <- reshape2::melt(to_plot[,c(1,4,6,5)])
to_plot$variable <- factor(to_plot$variable,levels=c("out","overlap"))
to_plot$position <- 100
to_plot$position[which(to_plot$variable=="overlap")] <- to_plot[which(to_plot$variable=="overlap"),"value"]
color <- setNames(c("red","gray"),c("overlap","out"))
ggplot(to_plot, aes(x = tissue, y = value, fill = variable)) +  
  geom_bar(stat = 'identity') +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank())+
  geom_text(data = to_plot,
            aes(label = sprintf("%.1f", value), y = position),
            color = "black", size = 5, vjust = 0.5)
