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
tissue <- "pancreas"
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
state_num <- 15
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum")) 
tissue_summary <- data.frame()
for(tissue in tissues){
  DML <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DML_delta01.txt"),header = T)
  DML <- DML[,c("chr","pos","fdr","diff")]
  DML$label <- paste0(DML$chr,"-",DML$pos)
  DML$end <- DML$pos+1
  DML <- as.data.table(DML)
  setDT(DML)
  setkey(DML,chr,pos,end)
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue,"_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_young$V2 <- chromHMM_young$V2+1
  bin <- chromHMM_young
  if(tissue %in% c("ovary","uterus","mammarygland")){
    bin <- bin[which(bin$V1%in% paste0("chr",c(1:19,"X"))),]
    bin$V1 <- factor(bin$V1,levels=paste0("chr",c(1:19,"X")))
  }else{
    bin$V1 <- factor(bin$V1,levels=paste0("chr",c(1:19,"X","Y")))
  }
  bin <- as.data.table(bin)
  setDT(bin)
  setkey(bin,V1,V2,V3)
  
  overlaps <- foverlaps(DML, bin, type = "any", nomatch = 0L)  
  
  increase <- overlaps[which(overlaps$fdr <0.05 & overlaps$diff>0),]
  decrease <- overlaps[which(overlaps$fdr <0.05 & overlaps$diff<0),]
  increase_table <- as.data.frame(table(increase$V4))
  decrease_table <- as.data.frame(table(decrease$V4))

  increase_table$condition <- "Hyper"
  decrease_table$condition <- "Hypo"

  
  increase_table$percentage <- increase_table$Freq/sum(increase_table$Freq)*100
  decrease_table$percentage <- decrease_table$Freq/sum(decrease_table$Freq)*100

  to_plot <- rbind(increase_table,decrease_table)
  to_plot$tissue <- tissue_label_change(tissue)
  to_plot$Var1 <- factor(to_plot$Var1,levels=paste0("E",1:15))
  # p <- ggplot(to_plot, aes(x = condition, y = percentage, fill = Var1)) +  
  #   geom_bar(stat = 'identity') +   
  #   theme_minimal() +   
  #   scale_fill_d3("category20") +
  #   theme(axis.title.x = element_blank(), 
  #         axis.text.x = element_text(angle = 45, hjust = 1),
  #         text = element_text(size = 20),legend.title = element_blank())
  tissue_summary <- rbind(tissue_summary,to_plot)
}
to_plot <- tissue_summary
to_plot$Var1 <- as.character(to_plot$Var1)
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,paste0("E",1:15))
dictionary <- list("E1"=1, "E2"=2, "E3"=3,
                   "E4"=4, "E5"=5, "E6"=7,
                   "E7"=8, "E8"=6, "E9"=9,
                   "E10"=10,"E11"=11,"E12"=15,
                   "E13"=12,"E14"=13,"E15"=14)
keys <- names(dictionary)
values <- unlist(dictionary)
to_plot$Var1 <- values[match(to_plot$Var1, keys)]
to_plot$Var1 <- paste0("E",to_plot$Var1)
to_plot$Var1 <- factor(to_plot$Var1,levels=paste0("E",c(1:15)))
tissue_order <- c("Mammary Gland","Cecum","Thymus","Uterus","iWAT","Stomach","Skin","Spleen",
                  "Muscle","Bone Marrow","Liver","Ileum","Testis","Cortex","Jejunum","Tongue",
                  "Hippocampus","Colon","Bladder","Aorta","Cerebellum","Lung","Heart","Kidney","BAT","Ovary","Pancreas")
to_plot$tissue <- factor(to_plot$tissue,levels=tissue_order)
ggplot(to_plot[which(to_plot$condition=="Hyper"),], aes(x = tissue, y = percentage, fill = Var1)) +  
  geom_bar(stat = 'identity') +   
  theme_bw() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Percentage")+
  ggtitle(NULL) 
