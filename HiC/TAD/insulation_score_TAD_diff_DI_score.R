rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(limma)
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
tissue <- "lung"
resolution <- "10000"

boundary_strength_caculate <- function(tissue,resolution){
  boundary <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD.csv"))
  boundary$start <- boundary$start+1
  boundary$label <- paste0(boundary$chr,"-",boundary$start,"-",boundary$end)
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  
  samples <- search_table$sample_name[which(search_table$tissue==tissue)]
  resolution_label <- paste0(as.numeric(resolution)/1000,"kb")
  DI_delta_summary <- data.frame()
  bin_num <- 10
  for(sample in samples){
    DI_delta <- data.frame()
    DI <- read.table(paste0("data/samples/HiC/",tissue,"/TAD/fanc_DI/",sample,"_",resolution_label,".directionality_2mb.bed"))
    DI <- DI[which(DI$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
    DI$label <- paste0(DI$V1,"-",DI$V2,"-",DI$V3)  
    for(chr in paste0("chr",c(1:19,"X","Y"))){
      t_DI <- DI[which(DI$V1 == chr),]
      t_boundary <- boundary[which(boundary$chr == chr),]
      for(i in c(1:nrow(t_boundary))){
        row_indices <- which(t_DI$label == t_boundary$label[i])
        t_DI_matrix <- t_DI[c((row_indices-bin_num):(row_indices+bin_num-1)),]
        if(any(is.na(t_DI_matrix$V5))){
          delta=NA
        }else{
          if( length((which(t_DI_matrix$V5 > 0))) == 0 | which(t_DI_matrix$V5 > 0)[1]==1){
            delta=abs(min(t_DI_matrix$V5) - max(t_DI_matrix$V5))
          }else{
            delta=abs(min(t_DI_matrix$V5[1:(which(t_DI_matrix$V5 > 0)[1]-1)]) - max(t_DI_matrix$V5[(which(t_DI_matrix$V5 > 0)[1]):(bin_num*2)]))
          }
        }
        t_DI_delta <- data.frame(label=t_boundary$label[i],delta=delta)
        DI_delta <- rbind(DI_delta,t_DI_delta)
      }
    }
    colnames(DI_delta)[2] <- sample
    if(nrow(DI_delta_summary)==0){
      DI_delta_summary <- DI_delta
    }else{
      DI_delta_summary <- merge(DI_delta_summary,DI_delta,by="label")
    }
  }
  rownames(DI_delta_summary) <- DI_delta_summary$label
  DI_delta_summary <- DI_delta_summary[,-1]
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  t_search_table$age[which(t_search_table$age=="3M")] <- "young"
  t_search_table$age[which(t_search_table$age=="24M")] <- "old"
  t_search_table$sample_name <- factor(t_search_table$sample_name, levels=colnames(DI_delta_summary))
  t_search_table <- t_search_table[order(t_search_table$sample_name),]
  t_search_table$age <- factor(t_search_table$age,levels=c("young","old"))
  
  DI_delta_summary <- na.omit(DI_delta_summary)
  DI_delta_summary <- DI_delta_summary[!apply(DI_delta_summary == 0, 1, all), ]
  
  design <- model.matrix(~age, t_search_table)
  fit <- lmFit(DI_delta_summary,design)
  fit <- eBayes(fit)
  out <- topTable(fit,coef=2,number = nrow(DI_delta_summary))
  out <- merge(DI_delta_summary,out,by = "row.names")
  out$Significant <- ifelse(out$adj.P.Val < 0.05 & abs(out$logFC) >= 0, 
                            ifelse(out$logFC > 0, "Up", "Down"), "Stable")
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(
    out, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (FDR)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," boundary strength"))+
    annotate("text", x = min(out$logFC), y = max(-log10(out$adj.P.Val)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$logFC), y = max(-log10(out$adj.P.Val)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  write.csv(out,paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD_diff_boundary_DI_score.csv"),row.names = F)
  return(p)
}
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus")
p_list <- list()
for(tissue in tissues){
  p_list[[tissue]] <- boundary_strength_caculate(tissue,resolution)
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_plot <- plot_a_list(p_list,no_of_rows = 3,no_of_cols = 4)
ggsave(paste0("result/HiC/all_tissues_insulation_score_",resolution,"_boundary_diff_DI_score.png"),combined_plot,width = 20,height = 15,type="cairo")
