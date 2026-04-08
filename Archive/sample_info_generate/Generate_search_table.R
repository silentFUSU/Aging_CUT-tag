rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
search_table <- data.frame(tissue = as.character(),
                           antibody = as.character(),
                           sample_name = as.character(),
                           mouse_ID = as.character(),
                           age = as.character())
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
    }
  }
  return(tissue_label)
}
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder","CB","jejunum","uterus","ovary")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody <- antibodys[j]
    if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
      df <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_10kb_bins.counts"),skip = 1)
    }else{
      df <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins.counts"),skip = 1)
    }
    t_sample_name <- colnames(df[,c(7:ncol(df))])
    pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
    t_sample_name <- gsub(pattern, "\\1", t_sample_name)
    t_mouse_ID <- c("101","110","102","111")
    t_tissue <- rep(tissue_label_change(tissue),4)
    t_age <- c("3m","24m","3m","24m")
    t_antibody <- rep(antibody,4)
    if(tissue == "skin"){
      t_mouse_ID <- c("From bing","From bing","From bing","From bing")
    }else if(antibody=="H3K9me3" & tissue=="testis"){
      t_mouse_ID <- c("103","106","102","111")
    }else if(antibody=="H3K9me3" & tissue=="liver"){
      t_mouse_ID <- c("101","110","103","106")
    }else if(tissue=="ovary"){
      t_mouse_ID <- c("207","217","204","214")
      t_age <- c("24m","3m","24m","3m")
    }else if(tissue=="uterus"){
      t_mouse_ID <- c("214","199","215","200")
    }
    t_search_table <- data.frame(tissue = t_tissue,
                                 antibody = t_antibody,
                                 sample_name = t_sample_name,
                                 mouse_ID = t_mouse_ID,
                                 age = t_age
                                 )
    search_table <- rbind(search_table, t_search_table)
    }
}
write.csv(search_table,"data/samples/all/CUTTag_search_table.csv",row.names = F)

search_table <- data.frame(tissue = as.character(),
                           antibody = as.character(),
                           sample_name = as.character(),
                           mouse_ID = as.character(),
                           age = as.character())
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder","CB","jejunum","uterus","ovary")

for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- read.delim(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_1kb_bins.counts"),skip = 1)
  t_sample_name <- colnames(df[,c(7:ncol(df))])
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  t_sample_name <- gsub(pattern, "\\1", t_sample_name)
  t_mouse_ID <- c("101","110","102","111")
  t_tissue <- rep(tissue_label_change(tissue),4)
  t_age <- c("3m","24m","3m","24m")
  t_antibody <- rep("ATAC",4)
  if(tissue == "brain"){
    t_age <- c("3m","3m","24m","24m")
    t_mouse_ID <- c("97","103","106","112")
  }else if(tissue == "liver"){
    t_mouse_ID <- c("101","110","103","106")
  }else if(tissue == "ovary"){
    t_age <- c("24m","3m","24m","3m")
    t_mouse_ID <- c("207","217","204","214")
  }else if(tissue == "uterus"){
    t_mouse_ID <- c("214","199","215","200")
  }
  t_search_table <- data.frame(tissue = t_tissue,
                               antibody = t_antibody,
                               sample_name = t_sample_name,
                               mouse_ID = t_mouse_ID,
                               age = t_age
  )
  search_table <- rbind(search_table, t_search_table)
}
write.csv(search_table,"data/samples/all/ATAC_search_table.csv",row.names = F)


search_table <- data.frame(tissue = as.character(),
                           antibody = as.character(),
                           sample_name = as.character(),
                           mouse_ID = as.character(),
                           age = as.character())
sample_info <- read.csv("data/samples/RNA/sample_tissue_info.csv")
tissues <- unique(sample_info$TissueName)
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- sample_info[which(sample_info$TissueName==tissue),]
  df <- df[order(df$SampleID),]
  t_sample_name <- df$SampleID
  t_mouse_ID <- c("101","110","102","111")
  t_tissue <- rep(tissue_label_change(tissue),nrow(df))
  t_age <- df$Age
  t_antibody <- rep("RNA",nrow(df))
  if(tissue == "jejunum"){
    t_mouse_ID <- c("102","111","101","110")
  }else if(tissue == "Cecum"){
    t_mouse_ID <- c("101","103","111","110","109")
  }else if(tissue == "Colon"){
    t_mouse_ID <- c("110","101","102","111","109")
  }else if(tissue == "Testis"){
    t_mouse_ID <- c("101","110","103","111")
  }else if(tissue == "Kidney"){
    t_mouse_ID <- c("103","110","102","111")
  }else if(tissue == "Pancreas"){
    t_mouse_ID <- c("110","102","111")
  }else if(tissue == "FC"){
    t_mouse_ID <- c("101","109","103","111")
  }else if(tissue %in% c("Hip","Tongue","Bladder","Heart","Stomach","CB","Thymus","Aorta")){
    t_mouse_ID <- c("100","106","103","109")
  }else if(tissue == "Uterus"){
    t_mouse_ID <- c("214","215","199","200")
  }else if(tissue == "Spleen"){
    t_mouse_ID <- c("96","97","106","109")
  }else if(tissue == "Skin"){
    t_mouse_ID <- rep("From bing",nrow(df))
  }
  t_search_table <- data.frame(tissue = t_tissue,
                               antibody = t_antibody,
                               sample_name = t_sample_name,
                               mouse_ID = t_mouse_ID,
                               age = t_age)
  search_table <- rbind(search_table, t_search_table)
}
write.csv(search_table,"data/samples/all/RNA_search_table.csv",row.names = F)
