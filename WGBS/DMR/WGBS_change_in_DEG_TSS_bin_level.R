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
library(DSS)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
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
search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
colnames(search_table)[3] <- "sample"
p_list <- list()
i <- 1
for(tissue in tissues){
  DEG <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene.csv"))
  WGBS <- read.csv(paste0("data/samples/WGBS/",tissue,"/compress2bin/TSS_region_1kb_all_depth.csv"))
  WGBS <- merge(WGBS,search_table[,-1],by="sample")
  WGBS$percent <- WGBS$percent*100
  WGBS$age <- factor(WGBS$age, levels =c("3M","24M"))
  increase_gene <- DEG[which(DEG$Significant=="Up"),]
  decrease_gene <- DEG[which(DEG$Significant=="Down"),]
  WGBS_increase_gene <- WGBS[which(WGBS$gene_name %in% increase_gene$X),]
  WGBS_decrease_gene <- WGBS[which(WGBS$gene_name %in% decrease_gene$X),]
  t <- t.test(WGBS_increase_gene$percent[which(WGBS_increase_gene$age=="24M")],WGBS_increase_gene$percent[which(WGBS_increase_gene$age=="3M")])
  p_list[[i]] <- ggplot(WGBS_increase_gene, aes(x = age, y = percent,fill=sample)) +  
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot() +  
    theme_minimal() +
    theme(text = element_text(size = 20)) +
    labs(title = paste0(tissue_label_change(tissue),"\nGene Expression Increase Region"), x = NULL, y = "CG%") +
    annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
             hjust = 1.1, vjust = -1.1, size = 5, colour = "red")+
    annotate("text", x = Inf, y = Inf, label = paste("old mean =",  round(t$estimate[[1]],2)  ),   
             hjust = 1.1, vjust = 1.1, size = 5, colour = "red") +
    annotate("text", x = -Inf, y = Inf, label = paste("young mean =",  round(t$estimate[[2]],2)  ),   
             hjust = 0, vjust = 1.1, size = 5, colour = "red")
  i <- i+1
  t <- t.test(WGBS_decrease_gene$percent[which(WGBS_decrease_gene$age=="24M")],WGBS_decrease_gene$percent[which(WGBS_decrease_gene$age=="3M")])
  p_list[[i]] <- ggplot(WGBS_decrease_gene, aes(x = age, y = percent,fill=sample)) +  
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot() +  
    theme_minimal() +
    theme(text = element_text(size = 20)) +
    labs(title = paste0(tissue_label_change(tissue),"\nGene Expression Decrease Region"), x = NULL, y = "CG%") +
    annotate("text", x = Inf, y = -Inf, label = paste("p-value =",  format(t$p.value, scientific = TRUE, digits = 3)  ),   
             hjust = 1.1, vjust = -1.1, size = 5, colour = "red")+
    annotate("text", x = Inf, y = Inf, label = paste("old mean =",  round(t$estimate[[1]],2)  ),   
             hjust = 1.1, vjust = 1.1, size = 5, colour = "red") +
    annotate("text", x = -Inf, y = Inf, label = paste("young mean =",  round(t$estimate[[2]],2)  ),   
             hjust = 0, vjust = 1.1, size = 5, colour = "red")
  i <- i+1
  
}
combined_plot <- plot_a_list(p_list,no_of_rows = length(p_list)/2,no_of_cols = 2)
ggsave("result/WGBS/all_tissues_change_in_DEG_TSS.png",combined_plot,width = 14,height = 6*length(p_list)/2, type="cairo")
