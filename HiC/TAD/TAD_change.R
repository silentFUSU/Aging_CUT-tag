rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggrepel)
library(ggplot2)
options(scipen = 999)  
tissue <- "cecum"
resolution <- "20000"
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
dommaincaller <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  young_samples <- search_table$sample_name[which(search_table$age=="3M")]
  old_samples <- search_table$sample_name[which(search_table$age=="24M")]
  samples_list <- list(young=young_samples,old=old_samples)
  ages <- c("young","old")
  tad_summary <- data.frame()
  for(age in ages){
    for(i in c(1:length(samples_list[[age]]))){
      sample <- samples_list[[age]][i]
      df <- read.table(paste0("data/samples/HiC/",tissue,"/TAD/domaincaller/",sample,".allValidPairs.",resolution,".bed"))
      df <- df[which(df$V1 %in% c(paste0("chr",c(1:19,"X","Y")))),]
      if(nrow(tad_summary)==0){
        tad_summary <- data.frame(sample=sample,tissue=tissue,age=age,counts=nrow(df))
      }else{
        tad_summary <- rbind(tad_summary,data.frame(sample=sample,tissue=tissue,age=age,counts=nrow(df)))
      }
    }
  }
  tad_summary$age <- factor(tad_summary$age,levels=c("young","old"))
  color <- read.table("data/samples/7_distinct_color.txt")
  color <- setNames(color$V1[c(1,3)],c("young","old"))
  ggplot(tad_summary,aes(x=age,y=counts,color = age))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
    geom_text(aes(label = sample), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
    scale_color_manual(values=color)+
    ggtitle(paste0(tissue_label_change(tissue)," TAD counts"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("TAD counts")
}

dommaincaller <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  young_samples <- search_table$sample_name[which(search_table$age=="3M")]
  old_samples <- search_table$sample_name[which(search_table$age=="24M")]
  samples_list <- list(young=young_samples,old=old_samples)
  ages <- c("young","old")
  tad_summary <- data.frame()
  for(age in ages){
    for(i in c(1:length(samples_list[[age]]))){
      sample <- samples_list[[age]][i]
      df <- read.delim(paste0("data/samples/HiC/",tissue,"/TAD/arrowhead/",sample,"_",resolution,"_SCALE_output/",resolution,"_blocks.bedpe"))
      df <- df[-1,]
      df <- df[which(df$X.chr1 %in% c(paste0("chr",c(1:19,"X","Y")))),]
      if(nrow(tad_summary)==0){
        tad_summary <- data.frame(sample=sample,tissue=tissue,age=age,counts=nrow(df))
      }else{
        tad_summary <- rbind(tad_summary,data.frame(sample=sample,tissue=tissue,age=age,counts=nrow(df)))
      }
    }
  }
  tad_summary$age <- factor(tad_summary$age,levels=c("young","old"))
  color <- read.table("data/samples/7_distinct_color.txt")
  color <- setNames(color$V1[c(1,3)],c("young","old"))
  ggplot(tad_summary,aes(x=age,y=counts,color = age))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
    geom_text(aes(label = sample), position = position_jitter(width = 0.2), vjust = -1, size = 3) +
    scale_color_manual(values=color)+
    ggtitle(paste0(tissue_label_change(tissue)," TAD counts"))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("TAD counts")
}

insulation_score <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  young_samples <- search_table$sample_name[which(search_table$age=="3M")]
  old_samples <- search_table$sample_name[which(search_table$age=="24M")]
  samples_list <- list(young=young_samples,old=old_samples)
  ages <- c("young","old")
  tad_summary <- data.frame()
  for(age in ages){
    for(i in c(1:length(samples_list[[age]]))){
      sample <- samples_list[[age]][i]
      df <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",sample,"_",resolution,"_tads.csv"))
      if(nrow(tad_summary)==0){
        tad_summary <- data.frame(sample=sample,tissue=tissue,age=age,counts=nrow(df))
      }else{
        tad_summary <- rbind(tad_summary,data.frame(sample=sample,tissue=tissue,age=age,counts=nrow(df)))
      }
    }
  }
  tad_summary$age <- factor(tad_summary$age,levels=c("young","old"))
  color <- read.table("data/samples/7_distinct_color.txt")
  color <- setNames(color$V1[c(1,3)],c("young","old"))
  p<- ggplot(tad_summary,aes(x=age,y=counts,color = age))+    
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7)+
    geom_text_repel(aes(label = sample), position = position_jitter(width = 0.2), size = 5) +
    # scale_color_manual(values=color)+
    ggtitle(paste0(tissue_label_change(tissue)," TAD counts"))+
    ylim(0,max(tad_summary$counts+1000))+
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+labs(fill = "", color = "") +ylab("Tad counts")
  return(p)
  print(p)
}
tissues <- sort(c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","skin","muscle","cecum","ileum"))
p_list <- list()
for(tissue in tissues){
  p_list[[tissue]] <- insulation_score(tissue,resolution)
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
combined_plot <- plot_a_list(p_list,no_of_rows = 4,no_of_cols = 4)
ggsave(paste0("result/HiC/all_tissues_insulation_score_",resolution,"_TAD_change.png"),combined_plot,width = 18,height = 24,type="cairo")
