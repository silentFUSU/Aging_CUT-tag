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
search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
search_table$tissue <- sapply(search_table$tissue,tissue_label_change)
search_table$age_label <- "Young"
search_table$age_label[which(search_table$age=="24m")] <- "Old"
search_table$rep <- "1"
search_table$rep[which(search_table$batch=="batch2")] <- "2"
search_table$rep[which(search_table$batch=="batch3")] <- "3"

search_table$age_short <- "3"
search_table$age_short[which(search_table$age=="24m")] <- "24"
search_table$sample_label <- paste0(search_table$tissue,"-",search_table$antibody,"-",search_table$age_label,search_table$rep)
df <- as.data.frame(table(search_table$sample_label))
# write.csv(search_table,"data/samples/sample_info/ATAC.csv")
search_table <- search_table[which(!search_table$tissue %in% c("Mef","Mef_ezh2_inhibit")),]
data_path <- read.table("data/samples/all/CUTTAG_ATAC_datapath.txt")
data_path <- data_path[which(data_path$V1 %in% search_table$sample_name),]
result <- merge(search_table, data_path, by.x = "sample_name", by.y="V1", sort = FALSE,all=T)
result <- result[order(result$tissue,result$sample_label),]
sample_data_path <- data.frame()
for(i in c(1:nrow(result))){
  t_result <- result[i,]
  if(!is.na(t_result$V2) & !is.na(t_result$V3)){
    value <- as.character(t_result[,c("V2","V3")])
    value <- sort(value)
    t_result$V2 <- value[1]
    t_result$V3 <- value[2]
    t_result$V2_label <- sub(".*/", "", t_result$V2)
    t_result$V3_label <- sub(".*/", "", t_result$V3)
    sample_data_path <- rbind(sample_data_path,t_result)
  }else{
    t_result$V2_label <- NA
    t_result$V3_label <- NA
    sample_data_path <- rbind(sample_data_path,t_result)
  }
}
write.csv(sample_data_path,"data/samples/sample_info/ATAC_data_path.csv",row.names = F)


