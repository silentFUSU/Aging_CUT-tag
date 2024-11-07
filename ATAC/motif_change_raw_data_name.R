rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
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
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
dir.create("data/samples/all/ATAC/motif_bg/all/")
dir.create("data/samples/all/ATAC/motif_bg/all/up")
dir.create("data/samples/all/ATAC/motif_bg/all/down")
for(tissue in tissues){
  for(condition in conditions){
    if(file.exists(paste0("data/samples/all/ATAC/motif_bg/",condition,"/",tissue,"/knownResults.txt"))){
      df <- read.delim(paste0("data/samples/all/ATAC/motif_bg/",condition,"/",tissue,"/knownResults.txt"))
      write.csv(df,paste0("data/samples/all/ATAC/motif_bg/all/",condition,"/",tissue_label_change(tissue),"_",condition,"_peaks_motif_enrichment.csv"),row.names = F)
    }
  }
}
