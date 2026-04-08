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
             "lung","muscle","skin","stomach","thymus","mammarygland")
tissue_summary <- data.frame()
for(tissue in tissues){
  annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
  split_names <- strsplit(annotation$label, "[:-]")
  annotation_df <- data.frame(
    chr = sapply(split_names, "[", 1),
    start = sapply(split_names, "[", 2),
    end = sapply(split_names, "[", 3),
    cluster = annotation$cluster
  )
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
  compartment <- compartment[,c("chr","start","end","young")]
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
  if(nrow(tissue_summary)==0){
    tissue_summary <- merged_counts
  }else{
    tissue_summary <- merge(tissue_summary,merged_counts,by="label",all=T)
  }
  }

to_plot <- tissue_summary
write.csv(to_plot, "data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table/kmeans_mean_compartmentB_percentage_detail.csv")
to_plot <- merge(to_plot,annotation,by="label")
to_plot <- to_plot[order(to_plot$cluster),]

rownames(to_plot) <- to_plot$label
to_plot <- to_plot[,-c(1,ncol(to_plot))]

rownames(annotation) <- annotation$label
annotation <- annotation[,-1,drop=F]
breaks <- c(seq(40, 100, length.out = 100))
color_palette <- colorRampPalette(c("white", "red"))(100)  
tissues_order <- c("Kidney","Muscle","Skin","Stomach","Heart","Hippocampus","Liver","Cortex","Cerebellum","Lung","Mammary Gland",
                   "Bone Marrow","Cecum","Colon","Thymus")
to_plot <- to_plot[,tissues_order]
annotation_color <- list(cluster=setNames(c("#f6416c", "#f8f3d4", "#ffde7d", "#00b8a9","grey"),c(paste0("kmeans",1:4),"Stable")))
pheatmap::pheatmap(to_plot,color = color_palette,cluster_cols = F,cluster_rows = F,show_rownames = F,annotation_row = annotation,annotation_colors = annotation_color)

annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
to_plot_box <- tissue_summary
to_plot_box <- merge(to_plot_box,annotation,by="label")
rownames(to_plot_box) <- to_plot_box$label
to_plot_box <- to_plot_box[,-1]
to_plot_box <- reshape2::melt(to_plot_box)
ggplot(to_plot_box, aes(x = cluster, y = value, fill = cluster)) +  
  geom_boxplot(outliers = F) + 
  theme_minimal() +   
  ylab("compartment B percentage")+
  ggtitle("Young")


