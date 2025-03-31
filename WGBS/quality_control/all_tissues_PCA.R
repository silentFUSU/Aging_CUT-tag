rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(Polychrome)
tissues <- c("mammarygland","lung","liver","kidney","ileum","Hip","skin","bonemarrow","jejunum","colon","ovary","CB","BAT","thymus","testis","heart","stomach","muscle","bladder","aorta","tongue","spleen","pancreas","brain","cecum","uterus","iWAT")
bin_size <- "1kb"

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

all_tissues_PCA <- function(tissues,bin_size){
  for(i in c(1:length(tissues))){
    tissue <- tissues[i]
    t_df <- read.csv(paste0("data/samples/WGBS/",tissue,"/compress2bin/",bin_size,"_bins_all_depth.csv"))
    t_df <- reshape2::dcast(t_df,label ~ sample, value.var = "percent")
    if(i == 1){
      df <- t_df
    }else{
      df <- merge(df,t_df,by="label")
    }
  }
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  rownames(df) <- df$label
  df <- df[,-1]
  df <- na.omit(df) 
  pca <- prcomp(t(df))
  saveRDS(pca,"data/samples/WGBS/1kb_bin_size_PCA.rds")
  to_plot <- data.frame(pca$x)
  to_plot$sample_name <- rownames(to_plot)
  to_plot <- merge(to_plot,search_table,by="sample_name")
  to_plot$age[which(to_plot$age == "3M")] <- "young"
  to_plot$age[which(to_plot$age == "24M")] <- "old"
  to_plot$age <- factor(to_plot$age,levels = c("young","old"))
  for(j in c(1:nrow(to_plot))){
    to_plot$tissue[j] <- tissue_label_change(to_plot$tissue[j])
  }
  
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  color_tissues <- sapply(sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
                                 "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")), tissue_label_change) 
  color <- read.table("data/samples/30_distinct_color.txt")
  color <- setNames(color$V1,color_tissues)
  
  p<- ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue, shape=age)) + 
    geom_point(size=5) +theme_bw()+
    scale_color_manual(values = color) +
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
    ggtitle(paste0("WGBS"))
    ggtitle(paste0("compress to ",bin_size," bin")) +
    geom_text_repel(  
      data = to_plot,  
      aes(x = PC1, y = PC2, label = mouse_ID, color = tissue),  
      size = 5,  
      box.padding = unit(0.35, "lines"),  
      point.padding = unit(0.3, "lines")  
    )
  return(p)
}
bin_size <- "1kb"
p <- all_tissues_PCA(tissues,bin_size)
p <- readRDS("data/samples/WGBS/1kb_bin_size_PCA_plot.rds")
ggsave("result/WGBS/all_tissue_PCA.png",p,width = 11,height = 8,type="cairo")
# saveRDS(p,"data/samples/WGBS/1kb_bin_size_PCA_plot.rds")

tissues <- c("lung","liver","ileum","kidney","Hip","bonemarrow","jejunum","colon","muscle","cecum")
per_tissue_PCA <- function(tissue,bin_size){
  df <- fread(paste0("data/samples/WGBS/",tissue,"/compress2bin/",bin_size,"_bins_all_depth.csv"))
  df <- data.table::dcast(df,label ~ sample, value.var = "percent")
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  setkey(df, label) 
  setattr(df, "row.names", df$label)  
  df[, label := NULL]  
  df <- df[complete.cases(df)]  
  df_matrix <- as.matrix(df) 
  pca <- prcomp(t(df_matrix))
  
  to_plot <- data.frame(pca$x)
  to_plot$sample_name <- rownames(to_plot)
  to_plot <- merge(to_plot,search_table,by="sample_name")
  to_plot$age[which(to_plot$age=="3M")] <- "Young"
  to_plot$age[which(to_plot$age=="24M")] <- "Old"
  to_plot$age <- factor(to_plot$age,levels=c("Young","Old"))
  to_plot$rownames <- paste0(to_plot$sample_name,"-",to_plot$mouse_ID,"-",to_plot$age)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  p <- ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
    geom_text_repel(  
      data = to_plot,  
      aes(x = PC1, y = PC2, label = rownames, color = age),  
      size = 5,  
      box.padding = unit(0.35, "lines"),  
      point.padding = unit(0.3, "lines")  
    ) +
    ggtitle(paste(tissue_label_change(tissue), "WGBS"))
  return(p)
}
p_list <- list()
for(tissue in tissues){
  p_list[[tissue]] <- per_tissue_PCA(tissue,bin_size)
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
combined_plot <- plot_a_list(p_list,2,5)
ggsave(paste0("result/WGBS/chrM_issue_tissues_PCA.png"),combined_plot,width = 24,height = 10,type="cairo")
