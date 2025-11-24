rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
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

tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")

tissue_summary <- list(up=data.frame(),down=data.frame())
for(tissue in tissues){
  sig_list <- list()
  search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  df <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3_diff_after_remove_batch_effect.csv"))
  sig_list[["up"]] <- df[which(df$Significant=="Up"),]
  sig_list[["down"]] <- df[which(df$Significant=="Down"),]
  for(condition in c("up","down")){
    if(nrow(sig_list[[condition]]) >= 20){
      tab <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3.counts"),header = T)
      tab_summary <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3.counts.summary"),header = T)
      rownames(tab) <- tab$Geneid
      counts = tab[,c(7:ncol(tab))]
      pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
      colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
      colnames(tab_summary) <- gsub(pattern,"\\1",colnames(tab_summary))
      
      t_search_table <- t_search_table[which(t_search_table$sample_name %in% colnames(counts)),]
      counts <- counts[,t_search_table$sample_name]
      tab_summary <- tab_summary[-2,t_search_table$sample_name]
      
      length_kb <- tab$Length / 1000  
      total_reads <- colSums(tab_summary)
      total_reads_million <- total_reads / 1e6  
      for (i in c(1:ncol(counts))) {  
        counts[[i]] <- (counts[[i]] / (length_kb * total_reads_million[i]))  
      }  
      
      RPKM <- counts
      young_cols <- RPKM[which(rownames(RPKM) %in% sig_list[[condition]]$Geneid), t_search_table$age=="3m"]
      young_cols$rowmeans <- rowMeans(young_cols)
      young_cols$label <- rownames(young_cols)
      young_cols$age <- "young"
      young_cols <- young_cols[,c("label","rowmeans","age")]
      
      old_cols <- RPKM[which(rownames(RPKM) %in% sig_list[[condition]]$Geneid), t_search_table$age=="24m"]
      old_cols$rowmeans <- rowMeans(old_cols)
      old_cols$label <- rownames(old_cols)
      old_cols$age <- "old"
      old_cols <- old_cols[,c("label","rowmeans","age")]
      
      t_tissue_summary <- rbind(young_cols,old_cols)
      t_tissue_summary$tissue <- tissue_label_change(tissue)
      
      tissue_summary[[condition]] <- rbind(tissue_summary[[condition]],t_tissue_summary)
    }
  }
}
#WGBS
WGBS_tissue_summary <- list(up=data.frame(),down=data.frame())
for(tissue in tissues){
  sig_list <- list()
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  df <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3_diff_after_remove_batch_effect.csv"))
  sig_list[["up"]] <- df[which(df$Significant=="Up"),]
  sig_list[["down"]] <- df[which(df$Significant=="Down"),]
  for(condition in c("up","down")){
    if(nrow(sig_list[[condition]]) >= 20){
      peaks <- as.data.table(sig_list[[condition]][,c(1:4)])
      setDT(peaks)
      setkey(peaks,Chr,Start,End)
      for(sample in t_search_table$sample_name){
        df_WGBS <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
        setDT(df_WGBS)
        setkey(df_WGBS,V1,V2,V3)  
        overlaps <- foverlaps(df_WGBS, peaks, type = "any", nomatch = 0L)  
        
        result <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(V5)), by = Geneid]
        result <- as.data.frame(result)
        result$methylation <- result$V4_sum/result$V5_sum
        result$age <- t_search_table$age[which(t_search_table$sample_name==sample)]
        result$sample <- sample
        result$tissue <- tissue_label_change(tissue)
        result <- result[,c("Geneid","methylation","age","sample","tissue")]
        WGBS_tissue_summary[[condition]] <- rbind(WGBS_tissue_summary[[condition]],result)
      }  
    }
  }
}
# saveRDS(WGBS_tissue_summary,"data/samples/ATAC/all/ATAC/WGBS_in_ATAC_diff_macs_young_old_narrowpeak_summits_spm3_diff_after_remove_batch_effect_tissue_peaks.rds")
WGBS_tissue_summary <- readRDS("data/samples/ATAC/all/ATAC/WGBS_in_ATAC_diff_macs_young_old_narrowpeak_summits_spm3_diff_after_remove_batch_effect_tissue_peaks.rds")
search_table <- read.csv("data/samples/all/WGBS_search_table.csv")

WGBS_mean_tissue_summary <- list()
for(condition in c("up","down")){
  summary_df <- WGBS_tissue_summary[[condition]] %>%
    group_by(Geneid,tissue, age) %>%
    summarise(value = mean(methylation, na.rm = TRUE))
  summary_df$age[which(summary_df$age=="3M")] <- "young"
  summary_df$age[which(summary_df$age=="24M")] <- "old"
  summary_df$antibody <- "WGBS"
  summary_df$value <- summary_df$value * 100
  WGBS_mean_tissue_summary[[condition]] <- summary_df
}

### histone
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
CUTTag_tissue_summary <- list(up=data.frame(),down=data.frame())
for(antibody in antibodys){
  for(tissue in tissues){
    sig_list <- list()
    search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
    t_search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
    df <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3_diff_after_remove_batch_effect.csv"))
    sig_list[["up"]] <- df[which(df$Significant=="Up"),]
    sig_list[["down"]] <- df[which(df$Significant=="Down"),]
    for(condition in c("up","down")){
      if(nrow(sig_list[[condition]]) >= 20){
        tab <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_in_ATAC_macs_young_old_narrowpeak_summits_spm3.counts"),header = T)
        tab_summary <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_in_ATAC_macs_young_old_narrowpeak_summits_spm3.counts.summary"),header = T)
        
        rownames(tab) <- tab$Geneid
        counts = tab[,c(7:ncol(tab))]
        pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
        colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
        colnames(tab_summary) <- gsub(pattern,"\\1",colnames(tab_summary))
        
        t_search_table <- t_search_table[which(t_search_table$sample_name %in% colnames(counts)),]
        counts <- counts[,t_search_table$sample_name]
        tab_summary <- tab_summary[-2,t_search_table$sample_name]
        
        length_kb <- tab$Length / 1000  
        total_reads <- colSums(tab_summary)
        total_reads_million <- total_reads / 1e6  
        for (i in c(1:ncol(counts))) {  
          counts[[i]] <- (counts[[i]] / (length_kb * total_reads_million[i]))  
        }  
        
        RPKM <- counts
        young_cols <- RPKM[which(rownames(RPKM) %in% sig_list[[condition]]$Geneid), t_search_table$age=="3m"]
        young_cols$rowmeans <- rowMeans(young_cols)
        old_cols <- RPKM[which(rownames(RPKM) %in% sig_list[[condition]]$Geneid), t_search_table$age=="24m"]
        old_cols$rowmeans <- rowMeans(old_cols)
        
        young_cols$label <- rownames(young_cols)
        young_cols$age <- "young"
        young_cols <- young_cols[,c("label","rowmeans","age")]
        
        old_cols$label <- rownames(old_cols)
        old_cols$age <- "old"
        old_cols <- old_cols[,c("label","rowmeans","age")]
        
        t_CUTTag_tissue_summary <- rbind(young_cols,old_cols)
        t_CUTTag_tissue_summary$tissue <- tissue_label_change(tissue)
        t_CUTTag_tissue_summary$antibody <- antibody
        
        CUTTag_tissue_summary[[condition]] <- rbind(CUTTag_tissue_summary[[condition]],t_CUTTag_tissue_summary)
      }
    }
  }
}
p_value_summary <- data.frame()
for(condition in c("up","down")){
  for(antibody in c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")){
    df <- CUTTag_tissue_summary[[condition]][which(CUTTag_tissue_summary[[condition]]$antibody==antibody),]
    young <- df[which(df$age=="young"),]
    old <- df[which(df$age=="old"),]
    test <- t.test(young$rowmeans,old$rowmeans)
    t_p_value_summary <- data.frame(condition=condition,antibody=antibody,p_value=test$p.value)
    p_value_summary <- rbind(p_value_summary,t_p_value_summary)
  }
}

for(condition in c("up","down")){
  df <- WGBS_tissue_summary[[condition]]
  young <- df[which(df$age=="3M"),]
  old <- df[which(df$age=="24M"),]

  test <- t.test(young$methylation,old$methylation)
  t_p_value_summary <- data.frame(condition=condition,antibody="WGBS",p_value=test$p.value)
  p_value_summary <- rbind(p_value_summary,t_p_value_summary)
}


for(condition in c("up","down")){
  tissue_summary[[condition]]$antibody <- "ATAC"
  WGBS_mean_tissue_summary[[condition]]<- WGBS_mean_tissue_summary[[condition]][,c("Geneid","value","age","tissue","antibody")]
  colnames(WGBS_mean_tissue_summary[[condition]]) <- c("label","rowmeans","age","tissue","antibody")
  
  to_plot <- rbind(tissue_summary[[condition]],WGBS_mean_tissue_summary[[condition]])
  to_plot <- rbind(to_plot,CUTTag_tissue_summary[[condition]])
  to_plot <- to_plot[which(to_plot$antibody %in% c("ATAC","H3K27ac","H3K4me1","H3K4me3")),]
  to_plot$age <- factor(to_plot$age,levels = c("young","old"))
  colnames(to_plot)[2] <- "value"
  p <- ggplot(to_plot, aes(x = antibody, y = value,fill=age)) +
    geom_boxplot(outlier.shape = NA) +
    theme_minimal()+  
    scale_fill_brewer(palette = "Pastel1") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
      axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
      axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
      axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
      legend.text = element_text(size = 12)
    )+ylim(0,15)
  ggsave(paste0("result/figures/ATAC_other_marks_change_boxplot_",condition,"_tissue_peak.pdf"),p,width = 8,height = 6)
}

for(condition in c("up","down")){
  to_plot <- rbind(tissue_summary[[condition]],WGBS_mean_tissue_summary[[condition]])
  to_plot <- rbind(to_plot,CUTTag_tissue_summary[[condition]])
  to_plot <- to_plot[which(to_plot$antibody %in% c("ATAC","H3K27ac","H3K4me1","H3K4me3")),]
  to_plot$age <- factor(to_plot$age,levels = c("young","old"))
  p <- ggplot(to_plot, aes(x = antibody, y = value,fill=age)) +
    geom_boxplot(outlier.shape = NA) +
    theme_minimal()+  
    scale_fill_brewer(palette = "Pastel1") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
      axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
      axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
      axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
      legend.text = element_text(size = 12)
    )+ylim(0,15)
  ggsave(paste0("result/figures/ATAC_other_marks_change_boxplot_",condition,".pdf"),p,width = 8,height = 6)
}

for(condition in c("up","down")){
  to_plot <- rbind(tissue_summary[[condition]],WGBS_mean_tissue_summary[[condition]])
  to_plot <- rbind(to_plot,CUTTag_tissue_summary[[condition]])
  to_plot <- to_plot[which(to_plot$antibody %in% c("WGBS")),]
  to_plot$age <- factor(to_plot$age,levels = c("young","old"))
  colnames(to_plot)[2] <- "value"
  p <- ggplot(to_plot, aes(x = antibody, y = value,fill=age)) +
    geom_boxplot(outlier.shape = NA) +
    theme_minimal()+  
    scale_fill_brewer(palette = "Pastel1") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
      axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
      axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
      axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
      legend.text = element_text(size = 12)
    )+ylim(0,100)
  ggsave(paste0("result/figures/ATAC_WGBS_change_boxplot_",condition,"_tissue_peak.pdf"),p,width = 3,height = 6)
}

