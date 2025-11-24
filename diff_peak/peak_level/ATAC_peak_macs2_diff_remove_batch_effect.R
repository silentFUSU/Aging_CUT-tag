rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
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
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
peak_preprocess_peak_level_remove_batch_effect <- function(tissue){
  search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
  tab = read.delim(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  counts <- counts[,search_table$sample_name]
  search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
  search_table$age <- factor(search_table$age, levels = c("3m","24m"))
  age <- as.character(search_table$age)
  batch <- as.character(search_table$batch)
  mouse_ID <- search_table$mouse_ID
  age[which(age=="3m")] <- "young"
  age[which(age=="24m")] <- "old"
  colnames(counts) <- paste0(colnames(counts),"-",age,"-",mouse_ID,"-",batch)
  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$year <- age
  y$samples$year <- factor(y$samples$year,c("young","old"))
  y$samples$batch <- search_table$batch
  y <- calcNormFactors(y)
  if(length(unique(y$samples$batch))==1){
    design <- model.matrix(~year, y$samples)
  }else{
    design <- model.matrix(~batch+year, y$samples)
  }
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = which(colnames(design) == "yearold"))
  tab<-tab[keep,]
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,
              "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC.old-young"=lrt$table$logFC)
  out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= log2(1.2), 
                            ifelse(out$`LogFC.old-young` > log2(1.2), "Up", "Down"), "Stable")
  write.csv(out,paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3_diff_after_remove_batch_effect.csv"),row.names = F)
  }
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  peak_preprocess_peak_level_remove_batch_effect(tissue)
}

### generate bed
dir.create("data/samples/ATAC/all/ATAC/macs2_diff_result/")
dir.create("data/samples/ATAC/all/ATAC/macs2_diff_result/stable/")
dir.create("data/samples/ATAC/all/ATAC/macs2_diff_result/up/")
dir.create("data/samples/ATAC/all/ATAC/macs2_diff_result/down/")
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3_diff_after_remove_batch_effect.csv"))
  increase <- df[which(df$LogFC.old.young >0 & df$FDR.old.young < 0.05),]
  decrease <- df[which(df$LogFC.old.young <0 & df$FDR.old.young < 0.05),]
  stable <- df[which(abs(df$LogFC.old.young) < 0.05 & df$FDR.old.young > 0.9),]
  if(nrow(increase) > 20){
    write.table(increase[,c("Chr","Start","End")],paste0("data/samples/ATAC/all/ATAC/macs2_diff_result/up/",tissue,"_summits_spm3.bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  }
  if(nrow(decrease) > 20){
    write.table(decrease[,c("Chr","Start","End")],paste0("data/samples/ATAC/all/ATAC/macs2_diff_result/down/",tissue,"_summits_spm3.bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  }
  write.table(stable[,c("Chr","Start","End")],paste0("data/samples/ATAC/all/ATAC/macs2_diff_result/stable/",tissue,"_summits_spm3.bed"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  }








