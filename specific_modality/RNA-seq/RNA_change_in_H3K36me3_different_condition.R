rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)

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
ref <- read.table("/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/TSS/refBed/mm10_refGene.bed")
ref$length <- ref$V3-ref$V2+1
ref <- ref %>%  
  group_by(V5) %>%  
  filter(length == max(length)) %>%  
  ungroup()  
ref <- as.data.table(ref)
setDT(ref)
setkey(ref, V1, V2, V3) 
tissue <- "mammarygland"

gene_expression_in_H3K36me3_different_condition <- function(tissue){
  tab <- read.delim(paste0("data/samples/RNA/",tissue,"/combined-chrM.counts"),skip=1)
  rownames(tab) <- tab$Geneid
  tab <- tab[,-1]
  colnames <- colnames(tab)[6:length(tab)]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+).*"
  colnames(tab)[6:length(tab)] <- gsub(pattern, "\\1", colnames(tab)[6:length(tab)] )
  counts <- tab[6:length(tab)]
  y= DGEList(counts=counts)
  RNA <- as.data.frame(cpm(y))
  
  histone <- read.csv(paste0("data/samples/",tissue,"/H3K36me3/H3K36me3_10kb_bins_diff.csv"))
  histone <- as.data.table(histone)
  setDT(histone)
  setkey(histone,Chr,Start,End)
  overlaps <- foverlaps(histone,ref, type = "any", nomatch = 0L) 
  result <- overlaps %>%  
    group_by(V5, Significant) %>% 
    summarise(count = n()) %>%  
    mutate(proportion = count / sum(count)) %>%  
    ungroup()
  
  result_cleaned <- result %>%  
    arrange(V5, desc(proportion), factor(Significant, levels = c("Down", "Up", "Stable"))) %>%  
    group_by(V5) %>%  
    slice_max(proportion, n = 1, with_ties = FALSE) %>%  
    ungroup() 
  
  increase <- result_cleaned$V5[which(result_cleaned$Significant=="Up")]
  decrease <- result_cleaned$V5[which(result_cleaned$Significant=="Down")]
  stable <- result_cleaned$V5[which(result_cleaned$Significant=="Stable")]
  
  RNA$H3K36me3_condition <- "Stable"
  RNA$H3K36me3_condition[which(rownames(RNA) %in% decrease)] <- "Down"
  RNA$H3K36me3_condition[which(rownames(RNA) %in% increase)] <- "Up"
  
  to_plot <- reshape2::melt(RNA)
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue_label_change(tissue)),]
  search_table$age <- factor(search_table$age, levels = c("3m","24m"))
  search_table <- search_table[order(search_table$age),]
  search_table$color <- c("#f38181", "#ff2e63", "#00adb5", "#3f72af")
  color <- setNames(search_table$color, search_table$sample_name)
  colnames(to_plot)[2] <- "sample_name"
  to_plot <- merge(to_plot,search_table,by="sample_name")
  to_plot$age <- factor(to_plot$age, levels=c("3m","24m"))
  to_plot$sample_name <- factor(to_plot$sample_name, levels = search_table$sample_name)
  to_plot$logcpm <- log2(to_plot$value+1)
  
  ggplot(to_plot, aes(x = H3K36me3_condition, y = logcpm,fill=sample_name)) +  
    geom_boxplot(outliers = F) +  
    scale_fill_manual(values = color) +
    labs(x = NULL,  
         y = "log2(CPM+1)") +
    theme_bw() +
    theme(  
      axis.text.x = element_text(angle = 45, hjust = 1)  
    ) +
    ggtitle(tissue_label_change(tissue),"Gene expression in H3K36me3 change regions")
  
  RNA$RNA_condition <- "Silent"
  diff_RNA <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene.csv"))
  RNA$RNA_condition[which(rownames(RNA) %in% diff_RNA$X[which(diff_RNA$Significant=="Up")])] <- "Up"
  RNA$RNA_condition[which(rownames(RNA) %in% diff_RNA$X[which(diff_RNA$Significant=="Stable")])] <- "Stable"
  RNA$RNA_condition[which(rownames(RNA) %in% diff_RNA$X[which(diff_RNA$Significant=="Down")])] <- "Down"
  
  H3K36me3_conditions <- c("Down","Stable","Up")
  to_plot <- data.frame()
  for(condition in  H3K36me3_conditions){
    t_to_plot <- RNA[which(RNA$H3K36me3_condition==condition),]  
    t_to_plot <- as.data.frame(table(t_to_plot$RNA_condition))
    t_to_plot$percent <- t_to_plot$Freq/sum(t_to_plot$Freq)*100
    t_to_plot$H3K36me3_condition <- condition
    to_plot <- rbind(to_plot,t_to_plot)
  }
  to_plot$H3K36me3_condition <- factor(to_plot$H3K36me3_condition,levels=c("Down","Stable","Up"))
  to_plot$Var1 <- factor(to_plot$Var1,levels=c("Silent","Stable","Up","Down"))
  ggplot(to_plot, aes(x = H3K36me3_condition, y = percent, fill = Var1)) +  
    geom_bar(stat = "identity") +  
    labs(x = "H3K36me3 condition", y = "Percentage", fill = "Gene expression condition") +  
    theme_minimal() +  
    theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle(paste0(tissue_label_change(tissue))) 
}
