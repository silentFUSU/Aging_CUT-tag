rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(readr) 
library(ggrepel)
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
search_table_CUTTAG <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
search_table_ATAC <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
search_table <- rbind(search_table_CUTTAG[,c("tissue","antibody","sample_name","age")],search_table_ATAC[,c("tissue","antibody","sample_name","age")])
antibodys <- c("ATAC","H3K27ac","H3K4me3")
qc_df <- data.frame()
for(i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  tsse <- read.table(paste0("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/tsse/",antibody,"/tsse.txt"), header = FALSE)
  pattern <- ".*(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  tsse$V2 <- gsub(pattern, "\\1", tsse$V2)
  colnames(tsse)[2] <- "sample_name"
  qc_df <- rbind(qc_df,tsse)
}

search_table <- search_table[which(search_table$antibody %in% antibodys),]
qc_df <- merge(search_table,qc_df[,c(2,3)],by="sample_name")
to_plot <- qc_df
to_plot$antibody <- factor(to_plot$antibody,levels=c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac","ATAC","RNA"))
# color <- read.table("data/samples/20_distinct_color.txt")
# color <- setNames(color$V1,c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac","ATAC","RNA"))
p <- ggplot(to_plot, aes(x = antibody, y = V3)) +
  geom_violin(fill="gray") +
  geom_boxplot(fill="white",width = 0.2,) +
  # scale_fill_manual(values = color) +
  labs(y = "TSSe") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
    axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
    axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
    axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
    legend.text = element_text(size = 12),  
    panel.background = element_blank(),  
    panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "lightgrey"),
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  ) + 
  ylim(0,20)
ggsave("result/Sup_figures/active_marks_TSSe.pdf",p,height=6,width = 5)

