rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
library(stringr)
library(dplyr)
library(ggrepel)
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


conditions <- c("up","down")
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
data_path <- "data/samples/ATAC/all/ATAC/snapatac2_macs/strict_stable_peaks_summits_spm3/"

motif_data_frame <- list(up=data.frame(),down=data.frame())
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(conditions))){
    condition <- conditions[j]
    if(file.exists(paste0(data_path,condition,"/enrichment_results_",tissue,"_summits_spm3.bed.csv"))){
      motif <- read.csv(paste0(data_path,condition,"/enrichment_results_",tissue,"_summits_spm3.bed.csv"))
      motif$adjusted.p.value[which(motif$log2.fold.change. < 0 )] <- 1
      motif$log2.fold.change.[which(motif$log2.fold.change. < 0 )] <- 0
      motif$adjusted.p.value[which(motif$fg_percent == 0 )] <- 1
      motif$log2.fold.change.[which(motif$fg_percent == 0 )] <- 0
      motif$log2.fold.change.[which(motif$fg_percent>0 & motif$bg_percent ==0)] <- max(motif$log2.fold.change.[is.finite(motif$log2.fold.change.)],na.rm = T)
      motif$id <- ifelse(
        grepl("\\(.*\\)", motif$id),
        sub(".*\\((.*?)\\).*", "\\1", motif$id), 
        sub("^M\\d+_2\\.00\\s*", "", motif$id) 
      )
      
      colnames(motif)[which(colnames(motif)=="adjusted.p.value")] <- tissue_label_change(tissue)
      if(i == 1){
        motif_data_frame[[condition]] <- motif[,c("id",tissue_label_change(tissue))]
      }else{
        motif_data_frame[[condition]] <- merge(motif_data_frame[[condition]],motif[,c("id",tissue_label_change(tissue))],by="id",all=T)
      }
    }else{
      motif <- data.frame(tissue=rep(1,790))
      colnames(motif) <- tissue_label_change(tissue)
      motif_data_frame[[condition]] <- cbind(motif_data_frame[[condition]],motif)
    }
  }
}
increase_count <- read.csv(paste0(data_path,"increase_common_motif_count.csv"))
decrease_count <- read.csv(paste0(data_path,"decrease_common_motif_count.csv"))
motifs <- unique(c(increase_count$id[which(increase_count$n >=6)],decrease_count$id[which(decrease_count$n>=6)]))

to_plot <- motif_data_frame[["up"]]
to_plot <- to_plot[which(to_plot$id %in% motifs),]
rownames(to_plot) <- to_plot$id
to_plot <- to_plot[,-1]

replace_zeros <- function(column) {
  non_zero_min <- min(column[column != 0], na.rm = TRUE)
  if (is.infinite(non_zero_min)) {
    non_zero_max <- 0
  }
  column[column == 0] <- non_zero_min
  return(column)
}
to_plot <- as.data.frame(apply(to_plot, 2, replace_zeros))
to_plot <- as.data.frame(apply(to_plot, 2, function(column) -log10(column)))

rowmeans <- as.data.frame(rowMeans(to_plot))
colnames(rowmeans) <- "means"
rowmeans <- rowmeans[order(rowmeans$means,decreasing = T),,drop=F]
to_plot <- to_plot[rownames(rowmeans),]

colmeans <- as.data.frame(colMeans(to_plot))
colnames(colmeans) <- "means"
colmeans  <- colmeans [order(colmeans$means,decreasing = T),,drop=F]
to_plot <- to_plot[,rownames(colmeans)]
breaks <- c(seq(-15,-(-log10(0.05)+0.0001), length.out = 80),seq(log10(0.05),-log10(0.05), length.out = 40),seq(-log10(0.05)+0.0001, 15, length.out = 80))
# breaks <- c(seq(-8,-0.51, length.out = 80),seq(-0.5,0.5, length.out = 40),seq(0.51, 8, length.out = 80))
color_palette <- c(
  colorRampPalette(c("blue","#defcf9"))(80),
  rep("white", 20), 
  rep("white", 20), 
  colorRampPalette(c("#ffe2e2","red"))(80) 
)
pheatmap::pheatmap(to_plot, cluster_rows =F,cluster_cols = F,color = color_palette,breaks=breaks,show_rownames = T,fontsize = 4)

to_plot <- motif_data_frame[["down"]]
rownames(to_plot) <- to_plot$id
to_plot <- to_plot[,-1]
replace_zeros <- function(column) {
  non_zero_min <- min(column[column != 0], na.rm = TRUE)
  if (is.infinite(non_zero_min)) {
    non_zero_max <- 0
  }
  column[column == 0] <- non_zero_min
  return(column)
}
to_plot <- as.data.frame(apply(to_plot, 2, replace_zeros))
to_plot <- as.data.frame(apply(to_plot, 2, function(column) -log10(column)))
to_plot <- to_plot[rownames(rowmeans),]

to_plot <- to_plot[,rownames(colmeans)]
to_plot <- -to_plot
pheatmap::pheatmap(to_plot, cluster_rows =F,cluster_cols = F,color = color_palette,breaks=breaks,show_rownames = T,fontsize = 4)

scatter_plot <- merge(increase_count[,c("id","n")],decrease_count[,c("id","n")],by="id",all=T)
scatter_plot[is.na(scatter_plot)] <- 0
colnames(scatter_plot)[2:3] <- c("increase_num","decrease_num")
scatter_plot$total_num <- scatter_plot$increase_num+scatter_plot$decrease_num

top_motif <- scatter_plot[which(scatter_plot$increase_num>=9 | scatter_plot$decrease_num>=9),] 
# top_motif <- scatter_plot[which(scatter_plot$total_num >=15),] 

p <- ggplot(scatter_plot, aes(x = increase_num, y = decrease_num, color=total_num)) +
  geom_point(size = 2) +
  labs(x = "# of tissues enriched, aging-up peaks", y = "# of tissues enriched, aging-down peaks",color = "Total tissue number") +
  theme_bw()+
  xlim(0,15)+
  ylim(0,15)+
  scale_color_gradient2(low = "#1fab89", mid = "#fff5a5", high = "red", midpoint = 8)+
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold")  
  )+
  geom_text_repel(data = top_motif, aes(x = increase_num, y = decrease_num, label = id), max.overlaps=100,
                  size = 5, 
                  nudge_y = 0.2)
ggsave("result/Sup_figures/ATAC_motif_scatter_plot.pdf",p,width = 7,height = 6)
