rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(ggsci)
library(gridExtra)
library(data.table)
tissue <- "ovary"
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
    }
  }
  return(tissue_label)
}
annotation_chromHMM_state <- function(tissue){
  DML <- readRDS(paste0("data/samples/WGBS/",tissue,"/DSS_table/dmlTest.sm.rds"))
  DML <- DML[,c("chr","pos","fdr","diff")]
  DML$label <- paste0(DML$chr,"-",DML$pos)
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  data_path <- paste0("data/samples/WGBS/",tissue,"/")
  for(sample in search_table$sample_name){
    cg <- fread(paste0(data_path,"/DSS_table/",sample,".txt"),sep="\t")
    cg$label <- paste0(cg$chr,"-",cg$pos)
    colnames(cg)[which(colnames(cg)=="N")] <- sample
    DML <- merge(DML,cg[,c(3,5)],by="label",all.x = TRUE)
  }
  saveRDS(DML,paste0("data/samples/WGBS/",tissue,"/DSS_table/dmlTest.sm.all.samples.rds"))
  rows_with_na <- apply(DML, 1, function(row) any(is.na(row)))  
  num_rows_with_na <- sum(rows_with_na)
  na_pie_plot <- data.frame(contidion=c("with_na","without_na"),count=c(num_rows_with_na,(nrow(DML)-num_rows_with_na)))
  na_pie_plot$percentage <- round(na_pie_plot$count / sum(na_pie_plot$count) * 100, 1) 
  na_pie_plot$label <- paste0(na_pie_plot$category, " (", na_pie_plot$percentage, "%)")  
  # ggplot(na_pie_plot, aes(x = "", y = count, fill = contidion)) +  
  #   geom_bar(stat = "identity", width = 1) +  
  #   coord_polar(theta = "y") +  
  #   theme_void() +  
  #   geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +  
  #   scale_fill_brewer(palette = "Set3")  
  
  DML_filtered <- na.omit(DML)
  file_dir <- paste0("result/all/ChromHMM/all_tissues/15_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_15_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  
  chromHMM_young_table <- as.data.frame(table(chromHMM_young$V4))
  chromHMM_young_table$Var1 <- factor(chromHMM_young_table$Var1,levels = paste0("E",c(1:15)))
  # ggplot(chromHMM_young_table, aes(x = "", y = Freq, fill = Var1)) +  
  #   geom_bar(stat = "identity", width = 1) +  
  #   coord_polar(theta = "y") +  
  #   theme_void() +  
  #   scale_fill_d3("category20")
  
  chromHMM_young <- data.table(chromHMM_young)
  DML_filtered <- data.table(DML_filtered)
  colnames(DML_filtered)[3] <- "start"
  DML_filtered$end <- DML_filtered$start
  setDT(chromHMM_young)  
  setDT(DML_filtered)  
  setkey(chromHMM_young, V1, V2, V3) 
  setkey(DML_filtered, chr, start, end) 
  overlaps <- foverlaps(DML_filtered, chromHMM_young, type = "any", nomatch = 0L)  
  increase <- overlaps[which(overlaps$fdr <0.05 & overlaps$diff>0),]
  decrease <- overlaps[which(overlaps$fdr <0.05 & overlaps$diff<0),]
  increase_table <- as.data.frame(table(increase$V4))
  decrease_table <- as.data.frame(table(decrease$V4))
  stable_table <- as.data.frame(table(overlaps$V4[which(overlaps$fdr>0.05)]))
  
  increase_table$condition <- "Hyper"
  decrease_table$condition <- "Hypo"
  stable_table$condition <- "Stable"
  
  increase_table$percentage <- increase_table$Freq/sum(increase_table$Freq)*100
  decrease_table$percentage <- decrease_table$Freq/sum(decrease_table$Freq)*100
  stable_table$percentage <- stable_table$Freq/sum(stable_table$Freq)*100
  
  to_plot <- rbind(increase_table,decrease_table)
  to_plot <- rbind(to_plot,stable_table)
  to_plot$Var1 <- factor(to_plot$Var1,levels=paste0("E",1:15))
  p <- ggplot(to_plot, aes(x = condition, y = percentage, fill = Var1)) +  
    geom_bar(stat = 'identity') +   
    theme_minimal() +   
    scale_fill_d3("category20") +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank())
  return(p)
}

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","colon","heart","Hip","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
p_list <- list()
for(i in c(1:length(tissues))){
  p_list[[i]] <- annotation_chromHMM_state(tissue = tissues[i])
}
saveRDS(p_list,"result/all/ChromHMM/until_ovary/15_until_ovary/p_list.rds")
# p_list <- readRDS("result/all/ChromHMM/until_ovary/15_until_ovary/p_list.rds")
# for(i in c(1:length(p_list))){
#   p_list[[i]] <- p_list[[i]] + ggtitle(tissue_label_change(tissues[i]))
# }
# combined_plot <- plot_a_list(p_list,no_of_rows = 3,no_of_cols = 5)
# ggsave("result/WGBS/all_tissue_DML_in_chromHMM_state.png",combined_plot,width = 15,height = 15,type="cairo")
