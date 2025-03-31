rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(edgeR)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
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
tissue <- "lung"
state_num <- 14
ATAC_in_chromHMM_state <- function(tissue,state_num){
  file_dir <- paste0("result/all/ChromHMM/all_tissues_previous/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  # chromHMM_young$V2 <- chromHMM_young$V2+1
  chromHMM_young <- data.table(chromHMM_young)
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_old[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_old <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_old <- chromHMM_old[which(chromHMM_old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  # chromHMM_old$V2 <- chromHMM_old$V2+1
  chromHMM_old <- data.table(chromHMM_old)
  
  search_table <- read.csv("data/samples/all/ATAC_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue_label_change(tissue)),]
  search_table$age <- factor(search_table$age, c("3m","24m"))
  search_table <- search_table[order(search_table$age),]
  
  tab <- read.delim(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_1kb_bins.counts"),skip=1)
  rownames(tab) <- tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(tab)[7:length(tab)] <-  gsub(pattern, "\\1", colnames(tab)[7:length(tab)])
  counts <- tab[7:length(tab)]
  y= DGEList(counts=counts)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  tab<-tab[keep,]
  tab$label <- paste0(tab$Chr,":",tab$Start,"-",tab$End)
  out <- cbind(tab[,"label",drop=F],cpm(y))
  p_list <- list()
  for(i in c(1:length(search_table$sample_name))){
    if(search_table$age[i] == "3m"){
      chromHMM <- chromHMM_young
    }else{
      chromHMM <- chromHMM_old
    }
    chromHMM$label <- paste0(chromHMM$V1,":",chromHMM$V2,"-",chromHMM$V3)
    chromHMM <- as.data.frame(chromHMM)
    to_plot <-merge(chromHMM[,c("label","V4")], out[,c("label",search_table$sample_name[i])],by="label")
    to_plot$V4 <- factor(to_plot$V4, levels = c(paste0("E",state_num:1)))
    colours <- read.table("data/samples/20_distinct_color.txt")
    colours <- setNames(colours$V1,c(paste0("E",1:state_num)))
    colnames(to_plot)[3] <- "CPM"
    to_plot$logCPM <- log2(to_plot$CPM)
    p_list[[i]] <- ggplot(to_plot, aes(x = V4, y = logCPM, fill=V4)) +  
      geom_boxplot() +
      scale_fill_manual(values = colours) +
      coord_flip() +
      theme_minimal()+
      theme(text = element_text(size = 20),legend.position = "none")+
      ggtitle(tissue_label_change(tissue),paste(search_table$sample_name[i],search_table$age[i]))+
      ylab("log2(CPM)") +
      xlab(NULL)
  }
  combined_plot <- plot_a_list(p_list,1,length(p_list))
  ggsave(paste0("result/",tissue,"/ATAC/ATAC_in_",state_num,"model_chromHMM.png"),combined_plot,height=8,width=6*length(p_list), type="cairo")
}

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")

for(tissue in tissues){
  ATAC_in_chromHMM_state(tissue,14)  
}
