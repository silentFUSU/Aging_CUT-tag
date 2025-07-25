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

data_path <-"data/samples/ATAC/ATAC_peak_from_LMJ/motif_snapatac2_cisbp/"
condition <-"down"
tissue <- "lung"
motif <- read.csv(paste0(data_path,"motif_all_peaks_bg/",condition,"/enrichment_results_",tissue,"_",toTitleCase(condition),"_sorted.bed.csv"))
motif <- motif[!grepl("\\)_\\(", motif$id), ]
motif$id <- sub("^M\\d+_2\\.00\\s*", "", motif$id)
motif <- motif[which(motif$adjusted.p.value <0.05 & motif$log2.fold.change. > 0 ),]
if(nrow(motif) >0 ){
  motif  <- motif[order(motif$adjusted.p.value),]
  motif <- motif[1:min(50,nrow(motif)),]
}

motif_PFMs <- read.csv(paste0(data_path,"motif_bg_db_from_PFM/",condition,"/enrichment_results_",tissue,"_",toTitleCase(condition),"_sorted.bed.csv"))
motif_PFMs <- motif_PFMs[!grepl("\\)_\\(", motif_PFMs$id), ]
motif_PFMs$id <- gsub(" Unknown","", motif_PFMs$id)
motif_PFMs <- motif_PFMs[which(motif_PFMs$adjusted.p.value <0.05 & motif_PFMs$log2.fold.change. > 0),]
if(nrow(motif_PFMs) >0 ){
  motif_PFMs <- motif_PFMs[order(motif_PFMs$adjusted.p.value),]
  motif_PFMs <- motif_PFMs[1:min(50,nrow(motif_PFMs)),]
}

motif <- as.character(motif$id)
motif_PFMs <- as.character(motif_PFMs$id)

overlap <- calculate.overlap(x = list(motif_PFMs = motif_PFMs, motif_meme = motif))
motif_PFMs <- overlap$a1
motif_meme <- overlap$a2
intersect <- overlap$a3
max_length <- max(length(motif_PFMs), length(motif_meme), length(intersect))

column1 <- c(motif_PFMs, rep("", max_length - length(motif_PFMs)))
column2 <- c(motif_meme, rep("", max_length - length(motif_meme)))
column3 <- c(intersect, rep("", max_length - length(intersect)))

result_table <- data.frame(
  `motif PFMs` = column1,
  `motif meme` = column2,
  `Intersect` = column3
)
tableGrob_obj <- tableGrob(result_table)
png(filename = paste0("result/all/ATAC/overlap_between_two_database/",condition,"/",tissue,"_",condition,"_motif.png"), width = 800, height = max_length*30)
grid.text(paste0(tissue_label_change(tissue)," ",condition), x = 0.5, y = unit(1, "npc") - unit(0.5, "lines"),
          gp = gpar(fontsize = 20, fontface = "bold"))
grid.draw(tableGrob_obj)
dev.off()



