rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
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
      tissue_label <- "Mammary Gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}

tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))

tissue <- "lung"
summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv"))
  rna <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  colnames(rna)[1] <- "Geneid"
  to_plot <- merge(rna[,c("Geneid","logFC","Significant")],df[,c("Geneid","LogFC.old.young","Significant")],by="Geneid")
  histone_increase_rna_increase <- nrow(to_plot[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Up"),])
  histone_increase_rna_decrease <- nrow(to_plot[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Up"),])
  histone_decrease_rna_increase <- nrow(to_plot[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Down"),])
  histone_decrease_rna_decrease <- nrow(to_plot[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Down"),])
  if(histone_increase_rna_increase > 5 & histone_increase_rna_decrease > 5 & histone_decrease_rna_increase > 5 & histone_decrease_rna_decrease > 5 & sum(histone_increase_rna_increase,histone_increase_rna_decrease,histone_decrease_rna_increase,histone_decrease_rna_decrease) > 40){
    chi_table <- matrix(c(histone_increase_rna_increase,histone_increase_rna_decrease,histone_decrease_rna_increase,histone_decrease_rna_decrease),nrow = 2, byrow = TRUE)
    dimnames(chi_table) <- list(
      RNA = c("increase", "decrease"),
      Histone = c("increase", "decrease")
    )
    chi_test <- chisq.test(chi_table)
    t_summary <- data.frame(tissue=tissue_label_change(tissue),chi_pvalue=chi_test$p.value)
  }else{
    t_summary <- data.frame(tissue=tissue_label_change(tissue),chi_pvalue=NA)
  }
  summary <- rbind(summary,t_summary)
}

summary$label <- "Insignificant"
summary$label[which(summary$chi_pvalue < 0.05)] <- "*"
summary$label[which(summary$chi_pvalue < 0.01)] <- "**"
summary$label[which(summary$chi_pvalue < 0.001)] <- "***"
chi_label_summary <- summary


tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))

tissue <- "lung"
summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_gene_TSS_diff_after_remove_batch_effect.csv"))
  rna <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  colnames(rna)[1] <- "Geneid"
  to_plot <- merge(rna[,c("Geneid","logFC","Significant")],df[,c("Geneid","LogFC.old.young","Significant")],by="Geneid")
  # histone_increase_rna_increase <- nrow(to_plot[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Up"),])
  histone_increase_rna_decrease <- nrow(to_plot[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Up"),])
  histone_decrease_rna_increase <- nrow(to_plot[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Down"),])
  # histone_decrease_rna_decrease <- nrow(to_plot[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Down"),])
  rna_decrease_histone_other <- nrow(rna[which(rna$Significant=="Down"),]) - histone_increase_rna_decrease
  rna_increase_histone_other <- nrow(rna[which(rna$Significant=="Up"),]) -  histone_decrease_rna_increase
  
  chi_table <- data.frame(histone_change=c(histone_decrease_rna_increase,histone_increase_rna_decrease),
                          histone_other=c(rna_increase_histone_other, rna_decrease_histone_other))
  rownames(chi_table) <- c("RNA_increase","RNA_decrease")
  chi_test <- chisq.test(chi_table)
  std_residuals <- chi_test$stdres
  t_summary <- data.frame(tissue=tissue_label_change(tissue),chi_pvalue=chi_test$p.value,residual=chi_test$stdres[1,1])
  # if(histone_decrease_rna_increase > 5 & histone_increase_rna_decrease > 5 &  rna_decrease_histone_other > 5 & rna_increase_histone_other & sum(histone_increase_rna_decrease,histone_decrease_rna_increase,rna_increase_histone_other,rna_decrease_histone_other) > 40){
  #   chi_table <- data.frame(histone_change=c(histone_decrease_rna_increase,histone_increase_rna_decrease),
  #                           histone_other=c(rna_increase_histone_other, rna_decrease_histone_other))
  #   rownames(chi_table) <- c("RNA_increase","RNA_decrease")
  #   chi_test <- chisq.test(chi_table)
  #   std_residuals <- chi_test$stdres
  #   t_summary <- data.frame(tissue=tissue_label_change(tissue),chi_pvalue=chi_test$p.value,residual=chi_test$stdres[1,1])
  # }else{
  #   fisher_table <- data.frame(histone_change=c(histone_decrease_rna_increase,histone_increase_rna_decrease),
  #                           histone_other=c(rna_increase_histone_other, rna_decrease_histone_other))
  #   rownames(fisher_table) <- c("RNA_increase","RNA_decrease")
  #   fisher_test <- fisher.test(fisher_table)
  #   std_residuals <- fisher_test$stdres
  #   t_summary <- data.frame(tissue=tissue_label_change(tissue),chi_pvalue=fisher_test$p.value,residual=NA)
  # }
  summary <- rbind(summary,t_summary)
}
summary$label <- ""
summary$label[which(summary$chi_pvalue < 0.05)] <- "*"
summary$label[which(summary$chi_pvalue < 0.01)] <- "**"
summary$label[which(summary$chi_pvalue < 0.001)] <- "***"
chi_label_summary <- summary
