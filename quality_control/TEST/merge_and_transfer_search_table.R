rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(stringr)
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
CUTTag <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff.csv")
CUTTag$antibody <- "CUT&Tag"
CUTTag <- CUTTag[,c(1,2,4,5)]
CUTTag <- unique(CUTTag) 
HiC <- read.csv("data/samples/all/HiC_search_table.csv")
WGBS <- read.csv("data/samples/all/WGBS_search_table.csv")
RNA <- read.csv("data/samples/all/RNA_search_table.csv")
df <- rbind(CUTTag,HiC[,c(1,2,4,5)],WGBS[,c(1,2,4,5)],RNA[,c(1,2,4,5)])

df$tissue <- sapply(df$tissue, tissue_label_change)  
df$label <- paste0(df$tissue,"-",df$antibody)
df$age[which(df$age=="3m" | df$age=="3M")] <- "Young"
df$age[which(df$age=="24m" | df$age=="24M")] <- "Old"
df$mouse_ID <- paste0(df$mouse_ID,"-",df$age)
label <- unique(df$label)
mouse <- unique(df$mouse_ID)
search_table <- data.frame(matrix(NA, nrow = length(mouse), ncol = length(label)))  
rownames(search_table) <- sort(mouse)  
colnames(search_table) <- label  
for(i in c(1:nrow(df))){
  search_table[df[i,"mouse_ID"],df[i,"label"]] <- 1
}
rowsum <- rowSums(search_table,na.rm = T)
search_table <- data.frame(RowSums = rowsum, search_table)

write.csv(search_table,"data/samples/all/all_samples_mouse_ID.csv",na = "")

#### get sample has problem
issue_mouse <- c("100","105","125","135","136","138","139","140","203","205","215","235","212","213","224","225","226","230","233","245","246","250")
CUTTag <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff.csv")
CUTTag <- CUTTag[which(CUTTag$mouse_ID %in% issue_mouse),]
write.csv(CUTTag,"data/samples/all/chrM_issue_CUTTag.csv")

WGBS <- read.csv("data/samples/all/WGBS_search_table.csv")
WGBS <- WGBS[which(WGBS$mouse_ID %in% issue_mouse),]
write.csv(WGBS,"data/samples/all/chrM_issue_WGBS.csv")


HiC <- read.csv("data/samples/all/HiC_search_table.csv")
HiC <- HiC[which(HiC$mouse_ID %in% issue_mouse),]
write.csv(HiC,"data/samples/all/chrM_issue_HiC.csv")

RNA <- read.csv("data/samples/all/RNA_search_table.csv")
RNA <- RNA[which(RNA$mouse_ID %in% issue_mouse),]
write.csv(RNA,"data/samples/all/chrM_issue_RNA.csv")
