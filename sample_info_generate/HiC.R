rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
tissue_label_change <- function(tissue){
  if(tissue=="brain"){
    tissue_label <- "Frontal Cortex"
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
search_table <- read.csv("data/samples/all/HiC_search_table.csv")
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
search_table <- search_table[which(search_table$tissue %in% tissues),]
search_table$tissue <- sapply(search_table$tissue,tissue_label_change)
search_table$age_label <- "Young"
search_table$age_label[which(search_table$age=="24M")] <- "Old"
search_table <- search_table[order(search_table$tissue,search_table$age),]
search_table$rep <- "1"
search_table$rep[which(search_table$batch=="batch2")] <- "2"
search_table$rep[which(search_table$batch=="batch3")] <- "3"

search_table$sample_label <- paste0(search_table$tissue,"-",search_table$antibody,"-",search_table$age_label,search_table$rep)
write.csv(search_table,"data/samples/sample_info/HiC.csv")
