rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(GenomeInfoDb)
library("GenomicRanges")
library(genomation)
library(data.table)
library(ggplot2)
library(stringr)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
cpg.file="~/ref_data/for_normal_mapping/mm10/cpgi.mm10.bed.txt"
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


cpgi_diff <- function(tissue){
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  cpg.shore.obj=readFeatureFlank(cpg.file,flank = 2000,feature.flank.name=c("CpGi","shores"))
  cpgi <- as.data.frame(cpg.shore.obj@listData$CpGi)
  cpgi$label <- paste(cpgi$seqnames,cpgi$start,cpgi$end,sep = "-")
  if(tissue %in% c("ovary","uterus","mammarygland")){
   cpgi <- cpgi[which(cpgi$seqnames %in% paste0("chr",c(1:19,"X"))),]
  }else{
    cpgi <- cpgi[which(cpgi$seqnames %in% paste0("chr",c(1:19,"X","Y"))),]
  }
  cpgi <- as.data.table(cpgi)
  setDT(cpgi)  
  setkey(cpgi,seqnames,start,end)
  summary <- data.frame()
  for(sample in t_search_table$sample_name){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    setDT(df)
    setkey(df,V1,V2,V3)  
    overlaps <- foverlaps(df, cpgi, type = "any", nomatch = 0L)  
    
    result <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(V5)), by = label]
    result <- as.data.frame(result)
    result$methylation <- result$V4_sum/result$V5_sum
    result <- result[,c("label","methylation")]
    colnames(result) <- c("label",sample)
    if(nrow(summary)==0){
      summary <- result
    }else{
      summary <- merge(summary,result,by="label")
    }
  }
  t_search_table$age[which(t_search_table$age=="3M")] <- "young"
  t_search_table$age[which(t_search_table$age=="24M")] <- "old"
  to_plot_box <- merge(to_plot_box,t_search_table,by.x="variable",by.y="sample_name")
  to_plot_box <- reshape2::melt(summary)
  to_plot_box$age <- factor(to_plot_box$age,levels = c("young","old"))
  ggplot(to_plot_box, aes(x = age, y = value, fill = age)) +  
    geom_boxplot(outliers = F) + 
    theme_minimal() +   
    ylab("DNA methylation")+
    ggtitle(tissue_label_change(tissue))
  
  }





