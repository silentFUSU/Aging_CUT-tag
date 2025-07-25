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
extract_before_bracket <- function(s) {  
  parts <- strsplit(s, "\\(")[[1]]  
  return(parts[1])  
}  
capitalize_first_lower_rest <- function(s) {
  if(nchar(s) > 0) { 
    paste0(toupper(substr(s, 1, 1)), tolower(substr(s, 2, nchar(s))))
  } else {
    s 
  }
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

data_path <-"data/samples/ATAC/ATAC_peak_from_LMJ/motif_homer_cisbp/motif_bg/"

tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
condition <-"down"
for(tissue in tissues){
  motif <- read.delim(paste0(data_path,condition,"/",tissue,"/knownResults.txt"))
  motif$Motif.Name <-  sapply(motif$Motif.Name, extract_before_bracket) 
  motif <- motif[which(motif$q.value..Benjamini. <0.05),]
  if(nrow(motif) >0 ){
    motif <- motif[order(motif$P.value),]
    motif <- motif[1:min(50,nrow(motif)),]
  }

  t_motif_snapatac2 <- read.csv(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/motif_snapatac2_cisbp/motif_bg_db_from_PFM/",condition,"/enrichment_results_",tissue,"_",toTitleCase(condition),"_sorted.bed.csv"))
  t_motif_snapatac2$id <- gsub(" Unknown", "", t_motif_snapatac2$id)
  t_motif_snapatac2 <- t_motif_snapatac2[which(t_motif_snapatac2$adjusted.p.value < 0.05 & t_motif_snapatac2$log2.fold.change. > 0),]
  if(nrow(t_motif_snapatac2) >0 ){
    t_motif_snapatac2 <- t_motif_snapatac2[order(t_motif_snapatac2$adjusted.p.value),]
    t_motif_snapatac2 <- t_motif_snapatac2[1:min(50,nrow(t_motif_snapatac2)),]
  }
  motif <- as.character(motif$Motif.Name)
  t_motif_snapatac2 <- as.character(t_motif_snapatac2$id)
  
  if(length(motif) > 0 | length(t_motif_snapatac2) > 0 ){
    overlap <- calculate.overlap(x = list(t_motif_snapatac2 = t_motif_snapatac2, motif = motif))
    only_motif_snapatac2 <- overlap$a1
    only_motif_homer <- overlap$a2
    intersect <- overlap$a3
    max_length <- max(length(only_motif_snapatac2), length(only_motif_homer), length(intersect))
    column1 <- c(only_motif_snapatac2, rep("", max_length - length(only_motif_snapatac2)))
    column2 <- c(only_motif_homer, rep("", max_length - length(only_motif_homer)))
    column3 <- c(intersect, rep("", max_length - length(intersect)))
    result_table <- data.frame(
      `motif snapatac2` = column1,
      `motif Homer` = column2,
      `Intersect` = column3
    )
    tableGrob_obj <- tableGrob(result_table)
    png(filename = paste0("result/all/ATAC/overlap_with_snapatac2_result/",condition,"/",tissue,"_",condition,"_motif.png"), width = 800, height = max_length*30)
    grid.text(paste0(tissue_label_change(tissue)," ",condition), x = 0.5, y = unit(1, "npc") - unit(0.5, "lines"),
              gp = gpar(fontsize = 20, fontface = "bold"))
    grid.draw(tableGrob_obj)
    dev.off()
  }
}

