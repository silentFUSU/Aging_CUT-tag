rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
library(tools)
library(kableExtra)  
library(VennDiagram) 
library(gridExtra)
library(grid)
library(png)

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

# increased_motif <- read.table("data/samples/ATAC/ATAC_peak_from_LMJ/motif_up_pvalue_matrix.txt")
# decreased_motif <- read.table("data/samples/ATAC/ATAC_peak_from_LMJ/motif_down_pvalue_matrix.txt")
data_path <-"data/samples/ATAC/ATAC_peak_from_LMJ/motif_snapatac2_cisbp/motif_all_peaks_bg/"

tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
condition <-"down"
for(tissue in tissues){
  motif <- read.csv(paste0(data_path,condition,"/enrichment_results_",tissue,"_",toTitleCase(condition),"_sorted.bed.csv"))
  motif <- motif[!grepl("\\)_\\(", motif$id), ]
  motif$id <- sub("^M\\d+_2\\.00\\s*", "", motif$id)
  motif <- motif[which(motif$adjusted.p.value <0.05 & motif$log2.fold.change. > 0),]
  if(nrow(motif) >0 ){
    motif  <- motif[order(motif$adjusted.p.value),]
    motif <- motif[1:min(50,nrow(motif)),]
  }
  
  t_motif_LMJ <- read.csv(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/from_LMJ/enrichment_",tissue,"_Peak_",toTitleCase(condition),"_",tissue,".csv"))
  t_motif_LMJ <- t_motif_LMJ[!grepl("\\)_\\(", t_motif_LMJ$id), ]
  t_motif_LMJ$id <- sub("^M\\d+_2\\.00\\s*", "", t_motif_LMJ$id)
  t_motif_LMJ <- t_motif_LMJ[which(t_motif_LMJ$adjusted.p.value < 0.05 & t_motif_LMJ$log2.fold.change. > 0),,drop=F]
  t_motif_LMJ <- t_motif_LMJ[order(t_motif_LMJ$adjusted.p.value),]
  t_motif_LMJ <- t_motif_LMJ[1:min(50,nrow(t_motif_LMJ)),,drop=F]
  
  
  t_motif_LMJ <- as.character(t_motif_LMJ$id)
  motif <- as.character(motif$id)
  
  if(length(motif) > 0 | length(t_motif_LMJ) > 0 ){
    overlap <- calculate.overlap(x = list(t_motif_LMJ = t_motif_LMJ, motif = motif))
    only_motif_LMJ <- overlap$a1
    only_motif_snapatac2_SZJ <- overlap$a2
    intersect <- overlap$a3
    max_length <- max(length(only_motif_LMJ), length(only_motif_snapatac2_SZJ), length(intersect))
    column1 <- c(only_motif_LMJ, rep("", max_length - length(only_motif_LMJ)))
    column2 <- c(only_motif_snapatac2_SZJ, rep("", max_length - length(only_motif_snapatac2_SZJ)))
    column3 <- c(intersect, rep("", max_length - length(intersect)))
    result_table <- data.frame(
      `motif LMJ` = column1,
      `motif snapatac2_SZJ` = column2,
      `Intersect` = column3
    )
    tableGrob_obj <- tableGrob(result_table)
    png(filename = paste0("result/all/ATAC/overlap_with_LMJ_result_snapatac2_all_peaks_bg/",condition,"/",tissue,"_",condition,"_motif.png"), width = 800, height = max_length*30)
    grid.text(paste0(tissue_label_change(tissue)," ",condition), x = 0.5, y = unit(1, "npc") - unit(0.5, "lines"),
              gp = gpar(fontsize = 20, fontface = "bold"))
    grid.draw(tableGrob_obj)
    dev.off()
  }
}



# for(tissue in tissues){
#   motif <- read.csv(paste0(data_path,condition,"/enrichment_results_",tissue,"_",toTitleCase(condition),"_sorted.bed.csv"))
#   motif <- motif[!grepl("\\)_\\(", motif$id), ]
#   motif$id <- sub("^M\\d+_2\\.00\\s*", "", motif$id)
#   motif <- motif[which(motif$adjusted.p.value <0.05 & motif$log2.fold.change. > 0),]
#   if(nrow(motif) >0 ){
#     motif  <- motif[order(motif$adjusted.p.value),]
#     motif <- motif[1:min(50,nrow(motif)),]
#   }
#   if(paste0("Peak_",toTitleCase(condition),"_",tissue) %in% colnames(motif_LMJ[[condition]])){
#     t_motif_LMJ <- motif_LMJ[[condition]][,paste0("Peak_",toTitleCase(condition),"_",tissue),drop=F]
#     t_motif_LMJ <- t_motif_LMJ[which(t_motif_LMJ[,1] > -log10(0.05)),,drop=F]
#     t_motif_LMJ <- t_motif_LMJ[1:min(50,nrow(t_motif_LMJ)),,drop=F]
#     t_motif_LMJ <- t_motif_LMJ[order(t_motif_LMJ[,1],decreasing = T),,drop=F]
#     t_motif_LMJ <- as.character(rownames(t_motif_LMJ))
#     
#     motif <- as.character(motif$id)
#     
#     if(length(motif) > 0 | length(t_motif_LMJ) > 0 ){
#       overlap <- calculate.overlap(x = list(t_motif_LMJ = t_motif_LMJ, motif = motif))
#       only_motif_LMJ <- overlap$a1
#       only_motif_snapatac2_SZJ <- overlap$a2
#       intersect <- overlap$a3
#       max_length <- max(length(only_motif_LMJ), length(only_motif_snapatac2_SZJ), length(intersect))
#       column1 <- c(only_motif_LMJ, rep("", max_length - length(only_motif_LMJ)))
#       column2 <- c(only_motif_snapatac2_SZJ, rep("", max_length - length(only_motif_snapatac2_SZJ)))
#       column3 <- c(intersect, rep("", max_length - length(intersect)))
#       result_table <- data.frame(
#         `motif LMJ` = column1,
#         `motif snapatac2_SZJ` = column2,
#         `Intersect` = column3
#       )
#       tableGrob_obj <- tableGrob(result_table)
#       png(filename = paste0("result/all/ATAC/overlap_with_LMJ_result_snapatac2_all_peaks_bg/",condition,"/",tissue,"_",condition,"_motif.png"), width = 800, height = max_length*30)
#       grid.text(paste0(tissue_label_change(tissue)," ",condition), x = 0.5, y = unit(1, "npc") - unit(0.5, "lines"),
#                 gp = gpar(fontsize = 20, fontface = "bold"))
#       grid.draw(tableGrob_obj)
#       dev.off()
#     }
#   }
# }
