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
    }
  }
  return(tissue_label)
}
search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
tissue <- "liver"
search_table <- search_table[which(search_table$tissue == tissue),]
df_list <- list()
for(i in c(1:nrow(search_table))){
  df_list[[i]] <- read.delim(paste0("data/raw_data/20240905_DYQ_005-018_WGBS/bed/",search_table$sample_name[i],"_CpG.bdg"),head=F)
  df_list[[i]]$depth <- df_list[[i]]$V5
  # df_list[[i]] <- df_list[[i]][which(df_list[[i]]$depth >15),]
  df_list[[i]]$percent <- df_list[[i]]$V4/df_list[[i]]$V5*100
  df_list[[i]]$label <- paste0(df_list[[i]]$V1,"-",df_list[[i]]$V2,"-",df_list[[i]]$V3)
  names(df_list)[i] <- search_table$sample_name[i]
}

young <- search_table$sample_name[which(search_table$age=="3M")]
young_list <- setNames(lapply(young, function(name) df_list[[name]]), young)
young_list <- lapply(names(young_list), function(name) {  
  young_list[[name]] %>%  
    filter(depth > 15) %>%  
    select(label, !!name := percent)
})  
young_df <- Reduce(function(x, y) merge(x, y, by = "label"), young_list)  

old <- search_table$sample_name[which(search_table$age=="24M")]
old_list <- setNames(lapply(old, function(name) df_list[[name]]), old)
old_list <- lapply(names(old_list), function(name) {  
  old_list[[name]] %>%  
    filter(depth > 15) %>%  
    select(label, !!name := percent)
})  
old_df <- Reduce(function(x, y) merge(x, y, by = "label"), old_list)  

par(cex.lab = 1.5, cex.axis = 1.2)  
smoothScatter(as.numeric(young_df[,3]) ~ as.numeric(young_df[,2]),xlab = colnames(young_df)[2],ylab = colnames(young_df)[3],main = "Young")
abline(a = 0, b = 1, col = "red", lty = 2)  
mtext(paste0(colnames(young_df)[2]," = ", round(mean(young_df[,2]),2),"%"), side = 3, line = -2, adj = 0.05,  cex = 1.2)  
mtext(paste0(colnames(young_df)[3]," = ", round(mean(young_df[,3]),2),"%"), side = 3, line = -3.5, adj = 0.05,  cex = 1.2)  
mtext(paste0("r = ", round(cor(young_df[,2],young_df[,3]),2)), side = 1, line = -1.5, adj = 0.9,  cex = 1.2)  

par(cex.lab = 1.5, cex.axis = 1.2)  
smoothScatter(old_df[,3] ~ old_df[,2],xlab = colnames(old_df)[2],ylab = colnames(old_df)[3],main = "Old")
abline(a = 0, b = 1, col = "red", lty = 2)  
mtext(paste0(colnames(old_df)[2]," = ", round(mean(old_df[,2]),2),"%"), side = 3, line = -2, adj = 0.05,  cex = 1.2)  
mtext(paste0(colnames(old_df)[3]," = ", round(mean(old_df[,3]),2),"%"), side = 3, line = -3.5, adj = 0.05,  cex = 1.2)  
mtext(paste0("r = ", round(cor(old_df[,2],old_df[,3]),2)), side = 1, line = -1.5, adj = 0.9,  cex = 1.2)  

merge_list <- lapply(names(df_list), function(name) {  
  df_list[[name]] %>%  
    filter(depth > 15) %>%  
    select(label, !!name := percent)
}) 
merge_list <- setNames(merge_list, names(df_list))
merge_df <- Reduce(function(x, y) merge(x, y, by = "label",all=TRUE), merge_list)  
merge_df_to_plot <- reshape2::melt(merge_df)
merge_df_to_plot <- merge_df_to_plot %>%  
  filter(!is.na(value))  
search_table$age <- factor(search_table$age, levels= c("3M","24M"))
search_table <- search_table[order(search_table$age),]
colnames(search_table)[3] <- "variable"
variable_order <- paste0(search_table$variable,"-",search_table$mouse_ID,"-",search_table$age)
merge_df_to_plot <- merge(merge_df_to_plot, search_table, by = "variable" )
merge_df_to_plot$variable_label <- paste0(merge_df_to_plot$variable, "-", merge_df_to_plot$mouse_ID, "-", merge_df_to_plot$age)
merge_df_to_plot$variable_label <- factor(merge_df_to_plot$variable_label, levels = variable_order)
ggplot(merge_df_to_plot, aes(x = variable_label, y = value, fill= age)) +  
  geom_violin(adjust = 2.5) +          
  scale_fill_brewer(palette = "Pastel1") +
  geom_boxplot(width = 0.1, color = "black", fill = "white", outlier.shape = NA) +  
  theme_minimal()+
  ggtitle(tissue_label_change(tissue))+
  theme(text = element_text(size = 20),legend.position = "none")+
  labs(x = NULL,y = "CpG%") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


merge_df_pca <- Reduce(function(x, y) merge(x, y, by = "label"), merge_list)  
pca <- prcomp(t(merge_df_pca[,-1]))
merge_df_pca_to_plot <- data.frame(pca$x)
merge_df_pca_to_plot$sample_name <- rownames(merge_df_pca_to_plot)
colnames(search_table)[3] <- "sample_name"
merge_df_pca_to_plot <- merge(merge_df_pca_to_plot, search_table, by="sample_name")
merge_df_pca_to_plot$label <- paste0(merge_df_pca_to_plot$sample_name,"-",merge_df_pca_to_plot$mouse_ID,"-",merge_df_pca_to_plot$age)
merge_df_pca_to_plot$age <- factor(merge_df_pca_to_plot$age,levels = c("3M","24M"))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(merge_df_pca_to_plot, aes(x=PC1, y=PC2, color=age)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(  
    data = merge_df_pca_to_plot,  
    aes(x = PC1, y = PC2, label = label, color = age),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  ) +
  ggtitle(paste(tissue_label_change(tissue), "WGBS"))
