rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)
library(ggsignif)
library(deepToolsDownstream)
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
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))

antibody <- "H3K9me3"
p_list <- list()
for(tissue in tissues){
  se <- importCount(paste0("result/all/H3K27me3_domain/matrix_median/",tissue,"_",antibody,"_change_in_H3K27me3_domain.mat.gz"))
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
  t_tissue_summary <- data.frame()
  for(i in c(1:length(se@assays@data))){
    df <- se@assays@data[[i]]
    colmean<- as.data.frame(colMeans(df))
    pattern <- "(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
    names <- gsub(pattern, "\\1",names(se@assays@data)[i])
    colmean$position <- 1:3000
    colmean$sample <- names
    colmean$age <- search_table$age[which(search_table$sample_name==names)]
    colnames(colmean)[1] <- "value"
    t_tissue_summary <- rbind(t_tissue_summary,colmean)
  }
  t_tissue_summary <- t_tissue_summary %>%
    group_by(age, position) %>% 
    summarize(mean_value = mean(value, na.rm = TRUE))
  t_tissue_summary$tissue <- tissue_label_change(tissue)
  to_plot <- t_tissue_summary %>%
    group_by(age,position) %>%
    summarize(
      mean_value = mean(mean_value),
      se = sd(mean_value) / sqrt(n())
    )
  to_plot$age <- factor(to_plot$age,c("3m","24m"))
  color <- setNames(c("#e64b35","#3c5488"),c("3m","24m"))
  p_list[[tissue_label_change(tissue)]] <- ggplot(to_plot, aes(position,color=age,y=mean_value)) + 
    geom_line(size=1.2) +xlim(500,2500)+
    ylab("RPKM") +
    scale_color_manual(values=color)+
    theme_bw() +
    theme(text = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  }
names <- names(p_list)
names <- sort(names)
p_list <- p_list[names]
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect",axes = "collect_x",axis_titles = "collect")
}
combined_plot <- plot_a_list(p_list,no_of_rows = 27,no_of_cols = 1)
ggsave(paste0("result/Sup_figures/",antibody,"_in_H3K27me3_domains.pdf"),combined_plot,width = 8,height = 20)
