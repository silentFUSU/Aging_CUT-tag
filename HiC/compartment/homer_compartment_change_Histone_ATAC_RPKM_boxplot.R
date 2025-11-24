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
      tissue_label <- "Mammary gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
tissues <- c("brain","CB", "kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip",
             "mammarygland", "stomach", "thymus","skin","muscle","cecum","ileum","pancreas","spleen")
tissue_summary <- data.frame()
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K4me3","H3K4me1","H3K27ac","ATAC")
for(antibody in antibodys){
  for(tissue in tissues){
    if(antibody == "ATAC"){
      search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
      tab <- read.table(paste0("data/samples/ATAC/",tissue,"/",antibody,"/",tissue,"_compartment_50000.counts"),header = T) 
      tab_summary <- read.table(paste0("data/samples/ATAC/",tissue,"/",antibody,"/",tissue,"_compartment_50000.counts.summary"),header=T)
    }else{
      search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
      tab <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_compartment_50000.counts"),header = T)
      tab_summary <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_compartment_50000.counts.summary"),header = T)
    }
    rownames(tab) <- tab$Geneid
    counts = tab[,c(7:ncol(tab))]
    pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HM[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
    colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
    colnames(tab_summary) <- gsub(pattern,"\\1",colnames(tab_summary))
    t_search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
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
    
    compartment <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_50000.csv"))
    compartment$label <- paste0(compartment$chr,":",compartment$start,"-",compartment$end)
    compartment <- compartment[,c("label","condition")]
    young_cols <- merge(young_cols,compartment,by.x="row.names",by.y="label")
    old_cols <- merge(old_cols,compartment,by.x="row.names",by.y="label")
    young_cols <- young_cols %>%
      group_by(condition) %>%
      summarise(mean_rowmeans = mean(rowmeans))
    young_cols$age <- "young"
    
    old_cols <- old_cols %>%
      group_by(condition) %>%
      summarise(mean_rowmeans = mean(rowmeans))
    old_cols$age <- "old"
    
    t_tissue_summary <- rbind(young_cols,old_cols)
    t_tissue_summary$tissue <- tissue_label_change(tissue)
    t_tissue_summary$antibody <- antibody
    tissue_summary <- rbind(tissue_summary,t_tissue_summary)
  }
}

p_value_summary <- data.frame()
for(condition in c("A-A","A-B","B-A","B-B")){
  for(antibody in antibodys){
    t_df <- as.data.frame(tissue_summary[which(tissue_summary$antibody==antibody & tissue_summary$condition==condition),])   
    t_df <- reshape2::dcast(t_df,tissue~age,value.var="mean_rowmeans")
    test <- t.test(t_df$young,t_df$old,paired = T)
    t_p_value_summary <- data.frame(antibody=antibody,condition=condition,p_value=test$p.value)  
    p_value_summary <- rbind(p_value_summary,t_p_value_summary)
  }
}

for(condition in c("A-A","A-B","B-A","B-B")){
  to_plot <- tissue_summary[which(tissue_summary$condition==condition),]
  to_plot$age <- factor(to_plot$age,levels=c("young","old"))
  to_plot$antibody <- factor(to_plot$antibody,levels=c("H3K9me3","H3K27me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3","ATAC"))
  color <- setNames(c("#f39b7f","#4dbbd5"),c("young","old"))
  p <- ggplot(to_plot, aes(x = antibody, y = mean_rowmeans,fill=age)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw()+  
    scale_fill_brewer(palette = "Pastel1") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
      axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
      axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
      axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
      legend.text = element_text(size = 12)
    ) +
    ylab("RPKM")+
    ggtitle(condition)
  ggsave(paste0("result/Sup_figures/compartment_",condition,"_histone_ATAC_change.pdf"),p,width = 8,height = 6)
  
}



