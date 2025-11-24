rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
library(stringr)
library(dplyr)
library(ggplot2)
extract_before_bracket <- function(s) {  
  parts <- strsplit(s, "\\(")[[1]]  
  return(parts[1])  
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
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
replace_zeros <- function(column) {
  non_zero_min <- min(column[column != 0], na.rm = TRUE)
  if (is.infinite(non_zero_min)) {
    non_zero_max <- 0
  }
  column[column == 0] <- non_zero_min
  return(column)
}
DMR_motif_count_summary <- data.frame()
DMR_motif_summary <- list(increase=data.frame(),decrease=data.frame())
DMR_motif_data_frame <- list(increase=data.frame(),decrease=data.frame())
conditions <- c("increase","decrease")
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
data_path <- "data/samples/WGBS/all/snapatac2/mutual_bg/"
ATAC_data_path <- "data//samples/ATAC/ATAC_peak_from_LMJ/motif/snapatac2/strict_stable_background/"
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  DMR_motif_summary <- data.frame()
  ATAC_motif_summary <- data.frame()
  for(j in c(1:length(conditions))){
    condition <- conditions[j]
    if(condition == "increase"){
      condition_ATAC <- "down"
    }else{
      condition_ATAC <- "up"
    }
    if(file.exists(paste0(data_path,condition,"/enrichment_results_",tissue,"_DMR_",condition,"_delta01.bed.csv"))){
      DMR_motif <- read.csv(paste0(data_path,condition,"/enrichment_results_",tissue,"_DMR_",condition,"_delta01.bed.csv"))
      DMR_motif$adjusted.p.value[which(DMR_motif$log2.fold.change. <= 0 )] <- 1
      DMR_motif$adjusted.p.value[which(DMR_motif$fg_percent==0)] <- 1
      
      DMR_motif$id <- ifelse(
        grepl("\\(.*\\)", DMR_motif$id),
        sub(".*\\((.*?)\\).*", "\\1", DMR_motif$id), 
        sub("^M\\d+_2\\.00\\s*", "", DMR_motif$id) 
      )
      DMR_motif <- DMR_motif[which(DMR_motif$log2.fold.change. > 0 & DMR_motif[,9] < 0.05),]
      colnames(DMR_motif)[which(colnames(DMR_motif)=="adjusted.p.value")] <- "DMR"
      DMR_motif <- DMR_motif[,c("id","DMR")]
      DMR_motif$DMR <- replace_zeros(DMR_motif$DMR)
      if(condition=="increase"){
        DMR_motif$DMR <- -log10(DMR_motif$DMR)
      }else{
        DMR_motif$DMR <- log10(DMR_motif$DMR)
      }
      DMR_motif_summary <- rbind(DMR_motif_summary,DMR_motif)
    }
    
    if(file.exists(paste0(ATAC_data_path,condition_ATAC,"/enrichment_results_",tissue,".bed.csv"))){
      ATAC_motif <- read.csv(paste0(ATAC_data_path,condition_ATAC,"/enrichment_results_",tissue,".bed.csv"))
      ATAC_motif$adjusted.p.value[which(ATAC_motif$log2.fold.change. < 0 )] <- 1
      ATAC_motif$adjusted.p.value[which(ATAC_motif$fg_percent == 0 )] <- 1
      ATAC_motif$id <- ifelse(
        grepl("\\(.*\\)", ATAC_motif$id),
        sub(".*\\((.*?)\\).*", "\\1", ATAC_motif$id), 
        sub("^M\\d+_2\\.00\\s*", "", ATAC_motif$id) 
      )
      ATAC_motif <- ATAC_motif[which(ATAC_motif$log2.fold.change. > 0 & ATAC_motif$adjusted.p.value < 0.01),]
      colnames(ATAC_motif)[which(colnames(ATAC_motif)=="adjusted.p.value")] <- "DAR"
      ATAC_motif <- ATAC_motif[,c("id","DAR")]
      ATAC_motif$DAR <- replace_zeros(ATAC_motif$DAR)
      if(condition_ATAC=="up"){
        ATAC_motif$DAR <- -log10(ATAC_motif$DAR)
      }else{
        ATAC_motif$DAR <- log10(ATAC_motif$DAR)
      }
      ATAC_motif_summary <- rbind(ATAC_motif_summary,ATAC_motif)
      ATAC_motif_summary <- ATAC_motif_summary %>%
        group_by(id) %>%
        slice_max(order_by = abs(DAR), n = 1)
    }
  }
  to_plot <- merge(DMR_motif_summary,ATAC_motif_summary,by="id",all=T)
  to_plot[is.na(to_plot)] <- 0
  rownames(to_plot) <- to_plot$id
  to_plot <- to_plot[,-1]
  to_plot$sum <- rowSums(to_plot)
  to_plot <- to_plot[order(to_plot$sum,decreasing = T),]
  breaks <- c(seq(-15,-(-log10(0.05)+0.0001), length.out = 80),seq(log10(0.05),-log10(0.05), length.out = 40),seq(-log10(0.05)+0.0001, 15, length.out = 80))
  # breaks <- c(seq(-8,-0.51, length.out = 80),seq(-0.5,0.5, length.out = 40),seq(0.51, 8, length.out = 80))
  color <- c(
    colorRampPalette(c("blue","#defcf9"))(80),
    rep("white", 20), 
    rep("white", 20), 
    colorRampPalette(c("#ffe2e2","red"))(80) 
  )
  pheatmap::pheatmap(to_plot[,-ncol(to_plot)],cluster_rows = F,cluster_cols = F,breaks = breaks,color = color,show_rownames = F, main = tissue_label_change(tissue))
}

