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
  df <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta01.txt"),header = T)
  for(condition in c("up","down")){
    if(condition == "up"){
      t_df <- df[which(df$areaStat > 0),]
    }else{
      t_df <- df[which(df$areaStat < 0 ),]
    }
    t_tissue_summary <- data.frame(tissue=rep(tissue_label_change(tissue),2),
                                 age=c("young","old"),
                                 value=c(median(t_df$meanMethy2)*100,median(t_df$meanMethy1)*100),
                                 antibody=c("WGBS","WGBS"))
    tissue_summary[[condition]] <- rbind(tissue_summary[[condition]],t_tissue_summary)
  }
}


### histone
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1","ATAC")
tissues_filtered <- c("Mammary Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen","Pancreas","Ovary","BAT","Kidney")
CUTTag_tissue_summary <- list(up=data.frame(),down=data.frame())
for(antibody in antibodys){
  for(tissue in tissues){
    for(condition in c("up","down")){
      if(condition == "up"){
        condition_label <- "increase"
      }else{
        condition_label <- "decrease"
      }
      if(antibody != "ATAC"){
        search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
        t_search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
        tab <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_DMR_",condition_label,"_delta01.counts"),header = T)
        tab_summary <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_DMR_",condition_label,"_delta01.counts.summary"),header = T)
      }else{
        search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
        t_search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
        tab <- read.table(paste0("data/samples/ATAC/",tissue,"/",antibody,"/",antibody,"_DMR_",condition_label,"_delta01.counts"),header = T)
        tab_summary <- read.table(paste0("data/samples/ATAC/",tissue,"/",antibody,"/",antibody,"_DMR_",condition_label,"_delta01.counts.summary"),header = T)
      }
      # print(paste0(tissue," ",condition," ",nrow(tab)))
      if(tissue %in% c("mammarygland","ovary","uterus")){
        tab <- tab[which(tab$Chr %in% paste0("chr",c(1:19,"X"))),]
      }
      # if (nrow(tab) < 1000) {
      #   next
      # }
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
      young_cols <- RPKM[, t_search_table$age=="3m"]
      young_cols$rowmeans <- rowMeans(young_cols)
      old_cols <- RPKM[, t_search_table$age=="24m"]
      old_cols$rowmeans <- rowMeans(old_cols)
      t_CUTTag_tissue_summary <- data.frame(tissue=rep(tissue_label_change(tissue),2),
                                            age=c("young","old"),
                                            value=c(median(young_cols$rowmeans),median(old_cols$rowmeans)),
                                            antibody=c(antibody,antibody))
      CUTTag_tissue_summary[[condition]] <- rbind(CUTTag_tissue_summary[[condition]],t_CUTTag_tissue_summary)
      }
    }
}
for(condition in c("up","down")){
  CUTTag_tissue_summary[[condition]] <- CUTTag_tissue_summary[[condition]][which(CUTTag_tissue_summary[[condition]]$tissue %in% tissues_filtered),]
}
p_value_summary <- data.frame()
for(condition in c("up","down")){
  for(antibody in c("ATAC","H3K27ac","H3K4me1","H3K4me3","H3K9me3","H3K27me3","H3K36me3")){
    df <- CUTTag_tissue_summary[[condition]][which(CUTTag_tissue_summary[[condition]]$antibody==antibody),]
    young <- df[which(df$age=="young"),]
    old <- df[which(df$age=="old"),]
    df <- merge(young,old,by="tissue")
    test <- t.test(df$value.x,df$value.y,paired=T)
    t_p_value_summary <- data.frame(condition=condition,antibody=antibody,p_value=test$p.value)
    p_value_summary <- rbind(p_value_summary,t_p_value_summary)
    }
}

for(condition in c("up","down")){
  to_plot <- rbind(tissue_summary[[condition]],CUTTag_tissue_summary[[condition]])
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
    ) +ylim(0,8)
  ggsave(paste0("result/figures/DMR_other_marks_change_boxplot_",condition,".pdf"),p,width = 4,height = 6)
}

### histone tissue-peak
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1","ATAC")
# tissues_filtered <- c("Mammary Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen","Pancreas","Ovary","BAT","Kidney")
CUTTag_tissue_summary <- list(up=data.frame(),down=data.frame())
for(antibody in antibodys){
  for(tissue in tissues){
    for(condition in c("up","down")){
      if(condition == "up"){
        condition_label <- "increase"
      }else{
        condition_label <- "decrease"
      }
      if(antibody != "ATAC"){
        search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
        t_search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
        tab <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_DMR_",condition_label,"_delta01.counts"),header = T)
        tab_summary <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_DMR_",condition_label,"_delta01.counts.summary"),header = T)
      }else{
        search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
        t_search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
        tab <- read.table(paste0("data/samples/ATAC/",tissue,"/",antibody,"/",antibody,"_DMR_",condition_label,"_delta01.counts"),header = T)
        tab_summary <- read.table(paste0("data/samples/ATAC/",tissue,"/",antibody,"/",antibody,"_DMR_",condition_label,"_delta01.counts.summary"),header = T)
      }
      # print(paste0(tissue," ",condition," ",nrow(tab)))
      if(tissue %in% c("mammarygland","ovary","uterus")){
        tab <- tab[which(tab$Chr %in% paste0("chr",c(1:19,"X"))),]
      }
      if (nrow(tab) < 1000) {
        next
      }
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
      young_cols <- RPKM[, t_search_table$age=="3m"]
      young_cols$rowmeans <- rowMeans(young_cols)
      old_cols <- RPKM[, t_search_table$age=="24m"]
      old_cols$rowmeans <- rowMeans(old_cols)
      
      young_cols$age <- "young"
      old_cols$age <- "old"
      young_cols$tissue <- tissue_label_change(tissue)
      old_cols$tissue <- tissue_label_change(tissue)
      young_cols$antibody <- antibody
      old_cols$antibody <- antibody
      young_cols <- young_cols[,c("rowmeans","tissue","age","antibody")]
      old_cols <- old_cols[,c("rowmeans","tissue","age","antibody")]
      
      t_CUTTag_tissue_summary <- rbind(young_cols,old_cols)
      CUTTag_tissue_summary[[condition]] <- rbind(CUTTag_tissue_summary[[condition]],t_CUTTag_tissue_summary)
    }
  }
}
p_value_summary <- data.frame()
for(condition in c("up","down")){
  for(antibody in c("ATAC","H3K27ac","H3K4me1","H3K4me3","H3K9me3","H3K27me3","H3K36me3")){
    df <- CUTTag_tissue_summary[[condition]][which(CUTTag_tissue_summary[[condition]]$antibody==antibody),]
    young <- df[which(df$age=="young"),]
    old <- df[which(df$age=="old"),]
    test <- t.test(young$rowmeans,old$rowmeans)
    t_p_value_summary <- data.frame(condition=condition,antibody=antibody,p_value=test$p.value)
    p_value_summary <- rbind(p_value_summary,t_p_value_summary)
  }
}

for(condition in c("up","down")){
  # to_plot <- rbind(tissue_summary[[condition]],CUTTag_tissue_summary[[condition]])
  to_plot <- CUTTag_tissue_summary[[condition]]
  to_plot <- to_plot[which(to_plot$antibody %in% c("ATAC","H3K27ac","H3K4me1","H3K4me3")),]
  to_plot$age <- factor(to_plot$age,levels = c("young","old"))
  colnames(to_plot)[1] <- "value"
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
    ) +ylim(0,8)
  ggsave(paste0("result/figures/DMR_other_marks_change_boxplot_",condition,"_tissue_peaks.pdf"),p,width = 4,height = 6)
}
