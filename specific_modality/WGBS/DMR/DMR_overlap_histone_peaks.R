rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(genomation)
library(methylKit)
library(ChIPseeker)
library(ggplot2)
library(stringr)
library(dplyr)

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}

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
DMR_annotation <- function(tissue,antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    peak <- paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W1000-G3000-E100.bed")
  }else{
    peak <- paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs_young_old_narrowpeak.bed")
  }
  p_list <- list()
  peak.shore.obj=readFeatureFlank(peak,flank=2000,feature.flank.name=c("peak region","near peak region"))
  DMR <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta0.txt"),header = T)
  colnames(DMR)[1:3] <- c("chr","start","end")
  increase <- DMR[which(DMR$diff.Methy > 0),c(1:3)] 
  decrease <- DMR[which(DMR$diff.Methy < 0),c(1:3)]
  increase_anno=annotateWithFeatureFlank(as(increase,"GRanges"),
                                        peak.shore.obj$`peak region`,peak.shore.obj$`near peak region`,
                                               feature.name="peak region",flank.name="near peak region")
  decrease_anno=annotateWithFeatureFlank(as(decrease,"GRanges"),
                                         peak.shore.obj$`peak region`,peak.shore.obj$`near peak region`,
                                         feature.name="peak region",flank.name="near peak region")
  increase_summary <- data.frame(location = c("Peak region","Other"),
                                 percent = c(increase_anno@precedence[[1]],
                                             100 - increase_anno@precedence[[1]]))
  decrease_summary <- data.frame(location = c("Peak region","Other"),
                                 percent = c(decrease_anno@precedence[[1]],
                                             100 - decrease_anno@precedence[[1]]))
  increase_summary$condition <- "Hyper"
  decrease_summary$condition <- "Hypo"
  summary <- rbind(increase_summary,decrease_summary)
  summary$location <- factor(summary$location, levels = c("Other","Peak region"))
  color <- setNames(c("#009980","#838B8B"),c("Peak region","Other"))
  p_list[[1]] <- ggplot( summary, aes(x = condition, y = percent, fill = location)) +  
    geom_bar(stat = 'identity',color="white") +   
    theme_minimal() +   
    scale_fill_manual(values = color) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Count Proportion")+
    ggtitle(paste0(tissue_label_change(tissue)))
  
  increase <- cbind(increase,increase_anno@members)
  increase$length <- increase$end-increase$start+1
  decrease <- cbind(decrease,decrease_anno@members)
  decrease$length <- decrease$end-decrease$start+1
  
  increase_summary <- data.frame(location = c("Peak region","Other"),
                                length = c(sum(increase$length[which(increase$`peak region`==1)]),sum(increase$length[which(increase$`peak region`!=1)])))
  decrease_summary <- data.frame(location = c("Peak region","Other"),
                                 length = c(sum(decrease$length[which(decrease$`peak region`==1)]),sum(decrease$length[which(decrease$`peak region`!=1)])))
  increase_summary$condition <- "Hyper"
  decrease_summary$condition <- "Hypo"
  increase_summary$percent <- increase_summary$length/sum(increase_summary$length)*100
  decrease_summary$percent <- decrease_summary$length/sum(decrease_summary$length)*100
  summary <- rbind(increase_summary,decrease_summary)
  summary$location <- factor(summary$location, levels = c("Other","Peak region"))
  p_list[[2]] <- ggplot( summary, aes(x = condition, y = percent, fill = location)) +  
    geom_bar(stat = 'identity',color="white") +   
    theme_minimal() +   
    scale_fill_manual(values = color) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Coverage Proportion")+
    ggtitle(paste0(tissue_label_change(tissue)))
  return(p_list)
}
tissues <-sort(c("liver","lung","kidney","ileum","Hip","mammarygland","skin","bonemarrow","jejunum","colon","ovary","CB","BAT","thymus","testis"))

antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K4me3","H3K4me1","H3K27ac")
for(antibody in antibodys){
  p_list <- list()
  for(i in c(1:length(tissues))){
    tissue <- tissues[i]
    p_list[[i]] <- DMR_annotation(tissue,antibody)
  }
  combine_plot <- plot_a_list(p_list,no_of_rows = 3,no_of_cols = 5) + patchwork::plot_annotation(title = paste0(antibody," peak region"),theme = theme(plot.title = element_text(size = 40,hjust = 0.5)))  
  ggsave(paste0("result/WGBS/all_tissues_DMR_delta0_annotation_",antibody,"_peak_region.png"),combine_plot,width = 15,height = 15,type="cairo")
  
}
