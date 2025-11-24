rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
library(ggnewscale)
options(scipen=0)
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
conditions <- c("up","down")
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
data_path <- "data/samples/ATAC/ATAC_peak_from_LMJ/motif_snapatac2_cisbp/motif_bg/"
if(data_path=="data/samples/ATAC/ATAC_peak_from_LMJ/motif_snapatac2_cisbp/motif_bg/"){
  label <- "(using stable peaks as background)"
}

motif_summary <- data.frame()
for(tissue in tissues){
  t_motif_summary <- data.frame()
  for(condition in conditions){
    if(file.exists(paste0(data_path,condition,"/enrichment_results_",tissue,"_",str_to_title(condition),"_sorted.bed.csv"))){
      motif <- read.csv(paste0(data_path,condition,"/enrichment_results_",tissue,"_",str_to_title(condition),"_sorted.bed.csv"))
      motif <- motif[which(motif$log2.fold.change. > 0 & motif$adjusted.p.value < 0.05),]
      motif$id <- ifelse(
        grepl("\\(.*\\)", motif$id),
        sub(".*\\((.*?)\\).*", "\\1", motif$id), 
        sub("^M\\d+_2\\.00\\s*", "", motif$id) 
      )
      motif <- motif[,c("id","log2.fold.change.")]
      if(condition == "down"){
        motif$`log2.fold.change.` <- -motif$`log2.fold.change.`
      }
      t_motif_summary <- rbind(t_motif_summary,motif)
    }
  }
  t_motif_summary <- t_motif_summary %>%
    group_by(id) %>%                             
    slice_max(order_by = abs(`log2.fold.change.`), n = 1) %>%   
    ungroup()     
  colnames(t_motif_summary)[2] <- tissue_label_change(tissue)
  if(nrow(motif_summary)==0){
    motif_summary <- t_motif_summary
  }else{
    motif_summary <- merge(motif_summary,t_motif_summary,by="id",all=T)
  }
}
motif_summary[is.na(motif_summary)] <- 0

to_plot <- motif_summary
rownames(to_plot) <- to_plot$id
to_plot <- to_plot[,-1]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
row_means <- rowMeans(to_plot, na.rm = TRUE)
to_plot_row_mean <- data.frame(row=rownames(to_plot),mean=row_means)
to_plot_row_mean <- to_plot_row_mean[order(to_plot_row_mean$mean,decreasing = T),]
breaks <- c(seq(-1, -0.29, length.out = 40), seq(-0.3, 0.3, length.out = 20), seq(0.31, 1, length.out = 40))
to_plot <- to_plot[to_plot_row_mean$row,]
pheatmap::pheatmap(to_plot,breaks = breaks,color = color_palette,cluster_rows = F,show_rownames = F)
