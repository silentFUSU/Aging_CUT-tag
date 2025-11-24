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

for(tissue in tissues){
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue_label_change(tissue)),]
  tab <- read.table(paste0("data/samples/RNA/",tissue,"/counts/",tissue,"_20000_redundant_tads.counts"),header = T)
  tab_summary <- read.table(paste0("data/samples/RNA/",tissue,"/counts/",tissue,"_20000_redundant_tads.counts.summary"),header = T)
  rownames(tab) <- tab$Geneid
  counts = tab[,c(7:ncol(tab))]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HM[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
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
  
  TAD <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_20000_TAD_diff_larger_250000.csv"))
  TAD <- TAD %>%
    separate(X, into = c("Chr", "Start", "End"), sep = "-", convert = TRUE)
  TAD$Start <- as.numeric(TAD$Start) -1
  TAD$label <- paste0(TAD$Chr,":",TAD$Start,"-",TAD$End)
  # colnames(TAD)[1] <- "label"
  colnames(TAD)[which(colnames(TAD)=="Significant")] <- "condition"
  TAD <- TAD[,c("label","condition")]
  young_cols <- merge(young_cols,TAD,by.x="row.names",by.y="label")
  old_cols <- merge(old_cols,TAD,by.x="row.names",by.y="label")
  
  young_cols <- young_cols[,c("condition","rowmeans")]
  old_cols <- old_cols[,c("condition","rowmeans")]
  young_cols$age <- "young"
  old_cols$age <- "old"
  
  t_tissue_summary <- rbind(young_cols,old_cols)
  t_tissue_summary$tissue <- tissue_label_change(tissue)
  t_tissue_summary$antibody <- "RNA"
  tissue_summary <- rbind(tissue_summary,t_tissue_summary)
}

to_plot <- tissue_summary
to_plot$age <- factor(to_plot$age,levels=c("young","old"))
to_plot <- to_plot[which(to_plot$condition %in% c("Up","Down")),]
to_plot$condition <- factor(to_plot$condition,levels = c("Up","Stable","Down"))
color <- setNames(c("#f39b7f","#4dbbd5"),c("young","old"))
p <- ggplot(to_plot, aes(x = condition, y = rowmeans,fill=age)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw()+  
  ylab("RPKM")+
  scale_fill_manual(values=color) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12,face = "bold", color = "black"),  
    axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
    axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
    axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
    legend.text = element_text(size = 12)
  ) +ylim(0,1.5)
# ggsave("result/Sup_figures/TAD_change_RNA_change.pdf",p,width = 6,height = 8)
p_value_summary <- data.frame()
for(condition in c("Up","Stable","Down")){
  t_df <- as.data.frame(tissue_summary[which(tissue_summary$antibody=="RNA" & tissue_summary$condition==condition),])   
  test <- t.test(t_df$rowmeans[which(t_df$age=="young")],t_df$rowmeans[which(t_df$age=="old")])
  t_p_value_summary <- data.frame(antibody="RNA",condition=condition,p_value=test$p.value)  
  p_value_summary <- rbind(p_value_summary,t_p_value_summary)
}
