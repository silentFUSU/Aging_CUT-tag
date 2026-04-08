rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(data.table)  
library(DSS)
library(edgeR)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}

tissues <- c("liver","lung","kidney","ileum","Hip","mammarygland")
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

bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    return("10kb")
  }else{
    return("1kb")
  }
}

search_table <-read.csv("data/samples/all/WGBS_search_table.csv")
colnames(search_table)[3] <- "sample"
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")

for(tissue in tissues){
  p_list <- list()
  i=1
  if(tissue == "ileum"){
    antibodys <- c("H3K27me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
  }else{
    antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
  }
  for(antibody in antibodys){
    diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff.csv"))
    if(antibody %in% c("H3K27ac","H3K4me3","H3K4me1")){
      peak <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size(antibody),"_in_young_old_merge_macs_narrowpeak.bed"))
    }else{
      peak <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size(antibody),"_in_young_old_merge-W1000-G3000-E100.bed"))
    }
    
    diff <- diff[which(diff$Geneid %in% peak$V4),]
    df <- read.csv(paste0("data/samples/WGBS/",tissue,"/compress2bin/",bin_size(antibody),"_bins_all_depth.csv"))
    df <- df[which(df$total_V5>15),]
    df$percent <- df$percent * 100
    df <- dcast(df, label ~ sample, value.var = "percent")  
    df <- na.omit(df)
    t_search_table <- search_table[which(search_table$sample %in% colnames(df)[-1]),]
    t_search_table$sample <- factor(t_search_table$sample, levels = colnames(df)[-1])
    t_search_table <- t_search_table[order(t_search_table$sample),]
    age <- t_search_table$age
    rownames(df) <- df$label
    df <- df[,-1]
    df_means <- apply(as.matrix(df),1,function(row){
      tapply(row,age,mean)
    })
    df_means <- t(df_means) 
    df_means <- as.data.frame(df_means)
    df_means <- df_means[,c("3M","24M")]
    rownames(df_means) <- rownames(df)
    df_means$logFC <- log2(df_means$`24M`/df_means$`3M`)
    df_means$Geneid <- rownames(df_means)
    to_plot <- merge(diff[,c(1,15)],df_means[,c(3,4)],by="Geneid")
    colnames(to_plot)[c(2,3)] <- c("histone","WGBS")
    p_list[[i]] <- ggplot()+
      geom_point(data=to_plot, mapping=aes(x=histone,y=WGBS),color = "grey",alpha=0.5) +  
      geom_point(data=to_plot[which(to_plot[,2]<0 & to_plot[,3]>0),], mapping=aes(histone,WGBS),color = "#ff9a00") +
      geom_point(data=to_plot[which(to_plot[,2]>0 & to_plot[,3]<0),], mapping=aes(histone,WGBS),color = "#48466d") +
      labs(x=antibody,
           y="WGBS") +
      theme_bw()+theme(text = element_text(size = 18))+
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggtitle(paste0(tissue_label_change(tissue)," WGBS and ",antibody))+
      annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]>0 & to_plot[,3]>0),])),x = Inf, y = Inf,hjust = 1.1, vjust = 1.2,colour="#00b8a9",size=5)+
      annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]<0 & to_plot[,3]>0),])),x = -Inf, y = Inf,hjust = -0.1, vjust = 1.2,colour="#ff9a00",size=5)+
      annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]<0 & to_plot[,3]<0),])),x = -Inf, y = -Inf,hjust = -0.1, vjust = -1.2,colour="#f6416c",size=5)+
      annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]>0 & to_plot[,3]<0),])),x = Inf, y = -Inf,hjust = 1.1, vjust = -1.2,colour="#48466d",size=5)
    i <- i+1
    # par(cex.lab = 1.5, cex.axis = 1.2)  
    # smoothScatter(as.numeric(to_plot[,3]) ~ as.numeric(to_plot[,2]),xlab = colnames(to_plot)[2],ylab = colnames(to_plot)[3],main = "Young")
    # abline(a = 0, b = 0, col = "red", lty = 2)  
    # abline(v = 0, col = "red", lty = 2) 
    
  }
  combined_plot <- plot_a_list(p_list,no_of_cols = ceiling(length(p_list)/2), no_of_rows = 2)
  ggsave(paste0("result/WGBS/",tissue,"/WGBS_change_in_histone_change/plot/all_antibodys_scatter_plot.png"),width = 7*length(p_list)/2,height = 14,type="cairo")
  }
