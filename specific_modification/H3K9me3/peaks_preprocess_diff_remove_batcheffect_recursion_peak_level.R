rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
window_size="5000"
gap_size="10000"
antibody <- "H3K9me3"
tissue <- "tongue"
peak_preprocess_peak_level_remove_batch_effect <- function(tissue,antibody){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_recursion.counts"),skip=1)
  }else{
    tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_young_old_narrowpeak.counts"),skip=1)
  }
  # tab <- tab[-which(tab$Chr %in% "chrY"),]
  # tab <- tab[which(tab$Length>=100000),]
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  counts <- counts[,search_table$sample_name]
  search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
  search_table$age <- factor(search_table$age, levels = c("3m","24m"))
  
  age <- as.character(search_table$age)
  batch <- as.character(search_table$batch)
  mouse_ID <- search_table$mouse_ID
  age[which(age=="3m")] <- "young"
  age[which(age=="24m")] <- "old"
  colnames(counts) <- paste0(colnames(counts),"-",age,"-",mouse_ID,"-",batch)
  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$year <- age
  y$samples$year <- factor(y$samples$year,c("young","old"))
  y$samples$batch <- search_table$batch
  
  y <- calcNormFactors(y)
  design <- model.matrix(~batch+year, y$samples)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = which(colnames(design) == "yearold"))
  tab<-tab[keep,]
  
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
              "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC.old-young"=lrt$table$logFC)
  
  out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= log2(1.2), 
                            ifelse(out$`LogFC.old-young` > log2(1.2), "Up", "Down"), "Stable")
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    write.csv(out,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_recursion_diff_after_remove_batch_effect.csv"),row.names = F)
  }else{
    write.csv(out,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_young_old_narrowpeak_diff_after_remove_batch_effect.csv"),row.names = F)
  }
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(
      out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
      geom_point(aes(color = Significant, size = Length),alpha=0.2) +
      scale_color_manual(values = colour) +
      scale_size_continuous(range = c(1, 5)) +
      geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
      labs(x="log2(fold change)",
           y="-log10 (FDR)") +
      theme_bw()+
      theme(text = element_text(size = 20))+
      ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
      annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
      annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  out <- out[which(out$Length > 200000),]
  p <- ggplot(
    out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
    geom_point(aes(color = Significant, size = Length),alpha=0.2) +
    scale_color_manual(values = colour) +
    scale_size_continuous(range = c(1, 5)) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (FDR)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  return(p)
}
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K4me3","H3K4me1","H3K27ac")
p_list <- list()
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody <- antibodys[j]
    p_list[[j]] <- peak_preprocess_peak_level_remove_batch_effect(tissue,antibody)
  }
  combined_plot <- plot_a_list(p_list,2,3)
  ggsave(paste0("result/",tissue,"/all_diff_peak_level_volcano_plot_remove_batch_effect.png"),combined_plot,width = 15,height = 10,type="cairo")
}

for (i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  for(j in c(1:length(tissues))){
    tissue <- tissues[j]
    p_list[[j]] <- peak_preprocess_peak_level_remove_batch_effect(tissue,antibody)
  }
  combined_plot <- plot_a_list(p_list,4,7)
  # ggsave(paste0("result/all/diff/",antibody,"/all_tissues_diff_-W",window_size,"-G",gap_size,"-E100_peaks_volcano_plot_remove_batch_effect.png"),combined_plot,width = 40,height = 18,type="cairo")
  ggsave(paste0("result/all/diff/",antibody,"/all_tissues_diff_-W",window_size,"-G",gap_size,"-E100_recursion_peaks_volcano_plot_remove_batch_effect.png"),combined_plot,width = 40,height = 18,type="cairo")
  # ggsave(paste0("result/",tissue,"/all_diff_peak_level_volcano_plot_remove_batch_effect.png"),combined_plot,width = 15,height = 10,type="cairo")
}

p_list <- list()
for (i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  for(j in c(1:length(tissues))){
    tissue <- tissues[j]
    colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
    out <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_diff_after_remove_batch_effect.csv"))
    out <- out[which(out$Length > 100000),]
    p_list[[tissue]] <- ggplot(
      out, aes(x = `LogFC.old.young`, y = -log10(`FDR.old.young`))) +
      geom_point(aes(color = Significant, size = Length),alpha=0.2) +
      scale_color_manual(values = colour) +
      scale_size_continuous(range = c(1, 5)) +
      geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
      labs(x="log2(fold change)",
           y="-log10 (p-value)") +
      theme_bw()+
      theme(text = element_text(size = 20))+
      ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
      annotate("text", x = min(out$`LogFC.old.young`), y = max(-log10(out$`FDR.old.young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
      annotate("text", x = max(out$`LogFC.old.young`), y = max(-log10(out$`FDR.old.young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
    }
  combined_plot <- plot_a_list(p_list,4,7)
  ggsave(paste0("result/all/diff/",antibody,"/all_tissues_diff_-W",window_size,"-G",gap_size,"-E100_100kb_peaks_volcano_plot_remove_batch_effect.png"),combined_plot,width = 40,height = 18,type="cairo")
}


MA_plot <- function(tissue,antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_recursion_diff_after_remove_batch_effect.csv"))
  }else{
    df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_young_old_narrowpeak_diff_after_remove_batch_effect.csv"))
  }
  
  df <- df[,c("Length","Geneid","logCPM","LogFC.old.young","Significant")]  
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  df <- df[which(df$Length > 200000),]
  p <- ggplot(
    df, aes(x = `logCPM`, y = `LogFC.old.young`)) +
    geom_point(aes(color = Significant, size = Length),alpha=0.2) +
    scale_color_manual(values = colour) +
    labs(x="Log2(CPM)",
         y="Log2(Fold Change)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
    annotate("text", x = max(df$logCPM), y = min(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Down"),]), vjust = 0, hjust = 1,colour="blue",size=5)+
    annotate("text", x = max(df$logCPM), y = max(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Up"),]), vjust = 1, hjust = 1,colour="red",size=5)
  # df<- df[which(df$Length > 100000),]
  # p <- ggplot(
  #   df, aes(x = `logCPM`, y = `LogFC.old.young`)) +
  #   geom_point(aes(color = Significant, size = Length),alpha=0.2) +
  #   scale_color_manual(values = colour) +
  #   labs(x="Log2(CPM)",
  #        y="Log2(Fold Change)") +
  #   theme_bw()+
  #   theme(text = element_text(size = 20))+
  #   ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
  #   annotate("text", x = max(df$logCPM), y = min(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Down"),]), vjust = 0, hjust = 1,colour="blue",size=5)+
  #   annotate("text", x = max(df$logCPM), y = max(df$LogFC.old.young), label = nrow(df[which(df$Significant=="Up"),]), vjust = 1, hjust = 1,colour="red",size=5)
  return(p)
}
p_list <- list()

for(tissue in tissues){
  p_list[[tissue]] <- MA_plot(tissue, "H3K9me3")  
}
combined_plot <- plot_a_list(p_list,4,7)
ggsave(paste0("result/all/diff/",antibody,"/all_tissues_diff_-W",window_size,"-G",gap_size,"-E100_peaks_MA_plot_remove_batch_effect.png"),combined_plot,width = 40,height = 18,type="cairo")


plist <- list()
for(tissue in tissues){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_bedtools_diff_after_remove_batch_effect.csv"))
  }else{
    df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_young_old_narrowpeak_diff_after_remove_batch_effect.csv"))
  }
  df <- df[which(df$Length > 100000),]
  df <- df[,c("Geneid","logCPM","LogFC.old.young","Significant")]
  # df <- df[which(df$Significant !="Stable"),]
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  plist[[tissue]] <- ggplot(df, aes(x = logCPM, y = LogFC.old.young)) +  
    geom_bin2d(bins = 50) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
    labs(  
      x = "Log2(CPM)",  
      y = "Log2(Fold Change)"
    ) + 
    scale_fill_continuous(type = "viridis", trans = "log2") + 
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody," bin in peaks"))+
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")
}
combined_plot <- plot_a_list(plist,4,7)
ggsave(paste0("result/all/diff/",antibody,"/all_tissues_scatter_plot_young_old_merge-W5000-G10000-E100_bedtools_100kb_filtered_peaks.png"),combined_plot,width = 35,height = 20,type="cairo")


summary <- data.frame()
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv"))  
  df <- df[which(df$Length > 200000),]
  df <- df[,c("Geneid","LogFC.old.young","FDR.old.young")]
  colnames(df) <- c("Geneid",paste0(tissue,".LogFC"),paste0(tissue,".FDR"))
  if(nrow(summary)==0){
    summary <- df
  }else{
    summary <- merge(summary,df,by="Geneid",all=T)
  }
  }
write.csv(summary,"data/samples/all/H3K9me3/recursion_peaks_diff_table/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv")




