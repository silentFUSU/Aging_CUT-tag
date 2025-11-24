rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggsignif)
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
tissues <- c("bonemarrow","brain","CB","cecum","colon","Hip","kidney","liver","heart",
             "lung","muscle","skin","stomach","thymus","mammarygland","ileum","pancreas","spleen")
tissue_summary <- list(young=data.frame(),old=data.frame())
for(tissue in tissues){
  annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
  split_names <- strsplit(annotation$label, "[:-]")
  annotation_df <- data.frame(
    chr = sapply(split_names, "[", 1),
    start = sapply(split_names, "[", 2),
    end = sapply(split_names, "[", 3),
    cluster = annotation$cluster
  )
  annotation_df$cluster[which(annotation_df$cluster=="kmeans1" & annotation_df$cluster=="chrY")] <- "kmeans-chrY"
  
  random_regions1 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed")
  random_regions1$cluster <- "random whole genome"
  colnames(random_regions1)[1:3] <- c("chr","start","end")
  annotation_df <- rbind(annotation_df,random_regions1)
  
  random_regions2<- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed")
  random_regions2$cluster <- "random out of peak"
  colnames(random_regions2)[1:3] <- c("chr","start","end")
  annotation_df <- rbind(annotation_df,random_regions2)
  
  if(tissue %in% c("mammarygland","ovary","uterus")){ 
    annotation_df <- annotation_df[which(annotation_df$chr %in% paste0("chr",c(1:19,"X"))),]
  }
  annotation_df$start <- as.numeric(annotation_df$start)
  annotation_df$end <- as.numeric(annotation_df$end)
  annotation_df$label <- paste0(annotation_df$chr,":",annotation_df$start,"-",annotation_df$end)
  annotation_df <- as.data.table(annotation_df)
  setDT(annotation_df)
  setkey(annotation_df,chr,start,end)
  compartment <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_50000.csv"))
  compartment <- compartment[,c("chr","start","end","young","old")]
  compartment <- as.data.table(compartment)
  setDT(compartment)
  setkey(compartment,chr,start,end)
  overlaps <- foverlaps(compartment, annotation_df, type = "any", nomatch = 0L)  
  
  count_B_by_label <- overlaps[young == "B", .N, by = label]
  total_count_by_label <- overlaps[, .N, by = label]
  merged_counts <- merge(count_B_by_label, total_count_by_label, by = "label", suffixes = c("_B", "_total"),all=T)
  merged_counts[, proportion_B := N_B / N_total *100]
  merged_counts <- as.data.frame(merged_counts)
  merged_counts$proportion_B[is.na(merged_counts$proportion_B)] <- 0
  colnames(merged_counts)[4] <- tissue_label_change(tissue)
  merged_counts <- merged_counts[,c(1,4)]
  if(nrow(tissue_summary[["young"]])==0){
    tissue_summary[["young"]] <- merged_counts
  }else{
    tissue_summary[["young"]] <- merge(tissue_summary[["young"]],merged_counts,by="label",all=T)
  }
  
  count_B_by_label <- overlaps[old == "B", .N, by = label]
  total_count_by_label <- overlaps[, .N, by = label]
  merged_counts <- merge(count_B_by_label, total_count_by_label, by = "label", suffixes = c("_B", "_total"),all=T)
  merged_counts[, proportion_B := N_B / N_total *100]
  merged_counts <- as.data.frame(merged_counts)
  merged_counts$proportion_B[is.na(merged_counts$proportion_B)] <- 0
  colnames(merged_counts)[4] <- tissue_label_change(tissue)
  merged_counts <- merged_counts[,c(1,4)]
  if(nrow(tissue_summary[["old"]])==0){
    tissue_summary[["old"]] <- merged_counts
  }else{
    tissue_summary[["old"]] <- merge(tissue_summary[["old"]],merged_counts,by="label",all=T)
  }
}

annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
random_regions1 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed")
random_regions1$cluster <- "random whole genome"
random_regions1$label <- paste0(random_regions1$V1,":",random_regions1$V2,"-",random_regions1$V3)
rownames(random_regions1) <- random_regions1$label
random_regions1 <- random_regions1[,c("label","cluster"),drop=F]
annotation <- rbind(annotation,random_regions1)
random_regions2<- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed")
random_regions2$cluster <- "random out of peak"
random_regions2$label <- paste0(random_regions2$V1,":",random_regions2$V2,"-",random_regions2$V3)
rownames(random_regions2) <- random_regions2$label
random_regions2 <- random_regions2[,c("label","cluster"),drop=F]
annotation <- rbind(annotation,random_regions2)

to_plot <- data.frame()
for(age in c("young","old")){
  df <- tissue_summary[[age]]
  df <- merge(df,annotation,by="label")
  rownames(df) <- df$label
  df <- df[,-1]
  df <- reshape2::melt(df)
  df$age <- age
  to_plot <- rbind(to_plot,df)
  }
to_plot$cluster <- factor(to_plot$cluster,levels=c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
to_plot$age <- factor(to_plot$age,levels=c("young","old"))
to_plot$value[is.na(to_plot$value)] <- 0
p <-ggplot(to_plot, aes(x = cluster, y = value, fill = age)) +  
  geom_boxplot(outliers = F) + 
  theme_minimal() +   
  theme_bw() +
  ggtitle("Young sample B compartment") +
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14,angle = 90, hjust = 1),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  ) 
p
ggsave("result/Sup_figures/H3K9me3_kmeans_compartment_box_plot_young_old.pdf",p,height = 4,width = 8)
