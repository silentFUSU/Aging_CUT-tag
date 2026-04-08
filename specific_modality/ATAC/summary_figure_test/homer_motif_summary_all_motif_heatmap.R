rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)

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
motif_data_frame <- list(up=data.frame(),down=data.frame())
conditions <- c("up","down")
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
data_path <- "data/samples/ATAC/ATAC_peak_from_LMJ/motif_homer_default_database/motif_bg/"
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(conditions))){
    condition <- conditions[j]
    if(file.exists(paste0(data_path,condition,"/",tissue,"/knownResults.txt"))){
      motif <- read.delim(paste0(data_path,condition,"/",tissue,"/knownResults.txt"))
      motif$Motif.Name <-  sapply(motif$Motif.Name, extract_before_bracket)  
      motif$Motif.Name <- paste0(motif$Motif.Name,"-",motif$Consensus)
      motif <- motif[!duplicated(motif$Motif.Name),]
      colnames(motif)[which(colnames(motif)=="q.value..Benjamini.")] <- tissue_label_change(tissue)
      if(i == 1){
        motif_data_frame[[condition]] <- motif[,c("Motif.Name",tissue_label_change(tissue))]
      }else{
        motif_data_frame[[condition]] <- merge(motif_data_frame[[condition]],motif[,c("Motif.Name",tissue_label_change(tissue))],by="Motif.Name",all=T)
      }
    }
  }
}
increase_count <- read.csv("data/samples/ATAC/ATAC_peak_from_LMJ/motif_homer_default_database/motif_bg/up/all_tissues_ATAC_peaks_increase_motif_count.csv")
decrease_count <- read.csv("data/samples/ATAC/ATAC_peak_from_LMJ/motif_homer_default_database/motif_bg/down/all_tissues_ATAC_peaks_decrease_motif_count.csv")
motifs <- unique(c(increase_count$Motif.Name[which(increase_count$n >=5)],decrease_count$Motif.Name[which(decrease_count$n>=5)]))
to_plot <- motif_data_frame[["up"]]
to_plot <- to_plot[which(to_plot$Motif.Name %in% motifs),]
rownames(to_plot) <- to_plot$Motif.Name
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
breaks <- c(seq(-3,-(-log10(0.05)+0.0001), length.out = 80),seq(log10(0.05),-log10(0.05), length.out = 40),seq(-log10(0.05)+0.0001, 3, length.out = 80))

color_palette <- c(
  colorRampPalette(c("blue","#defcf9"))(80),
  rep("white", 20), 
  rep("white", 20), 
  colorRampPalette(c("#ffe2e2","red"))(80) 
)
pheatmap::pheatmap(to_plot, cluster_rows =F,cluster_cols = F,color = color_palette,breaks=breaks,show_rownames = T,fontsize = 8)

to_plot <- motif_data_frame[["down"]]
rownames(to_plot) <- to_plot$Motif.Name
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
pheatmap::pheatmap(to_plot, cluster_rows =F,cluster_cols = F,color = color_palette,breaks=breaks,show_rownames = T,fontsize = 8)


