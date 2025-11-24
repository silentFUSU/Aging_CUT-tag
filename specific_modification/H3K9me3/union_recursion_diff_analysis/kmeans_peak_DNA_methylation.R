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


tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
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
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue & search_table$age=="3M"),]
  summary <- data.frame()
  for(sample in t_search_table$sample_name){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    setDT(df)
    setkey(df,V1,V2,V3)  
    overlaps <- foverlaps(df, annotation_df, type = "any", nomatch = 0L)  
    
    result <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(V5)), by = label]
    result <- as.data.frame(result)
    result$methylation <- result$V4_sum/result$V5_sum * 100
    result <- result[,c("label","methylation")]
    colnames(result) <- c("label",sample)
    if(nrow(summary)==0){
      summary <- result
    }else{
      summary <- merge(summary,result,by="label")
    }
  }
  summary$mean <- rowMeans(summary[,-1])
  summary <- summary[,c("label","mean")]
  colnames(summary)[2] <- tissue_label_change(tissue)
  if(nrow(tissue_summary)==0){
    tissue_summary <- summary
  }else{
    tissue_summary <- merge(tissue_summary,summary,by="label",all=T)
  }
}

to_plot <- tissue_summary

to_plot <- merge(to_plot,annotation,by="label")
to_plot <- to_plot[order(to_plot$cluster),]

rownames(to_plot) <- to_plot$label
to_plot <- to_plot[,-c(1,ncol(to_plot))]

rownames(annotation) <- annotation$label
annotation <- annotation[,-1,drop=F]
breaks <- c(seq(40, 100, length.out = 100))
color_palette <- colorRampPalette(c("white", "red"))(100)  
tissues_order <- c("Kidney","Muscle","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung","Mammary Gland",
                   "Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
to_plot <- to_plot[,tissues_order]
annotation_color <- list(cluster=setNames(c("#f6416c", "#f8f3d4", "#ffde7d", "#00b8a9","grey"),c(paste0("kmeans",1:4),"Stable")))
pheatmap::pheatmap(to_plot,cluster_cols = F,cluster_rows = F,show_rownames = F,annotation_row = annotation,annotation_colors = annotation_color)

## box plot
annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
to_plot_box <- tissue_summary
to_plot_box <- merge(to_plot_box,annotation,by="label")
rownames(to_plot_box) <- to_plot_box$label
to_plot_box <- to_plot_box[,-1]
to_plot_box <- reshape2::melt(to_plot_box)
ggplot(to_plot_box, aes(x = cluster, y = value, fill = cluster)) +  
  geom_boxplot(outliers = F) + 
  theme_minimal() +   
  ylab("DNA methylation")+
  ggtitle("Young")

tissue_order <- c("Kidney","Muscle","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung",
                  "Mammary Gland","Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
to_plot_box$variable <- factor(to_plot_box$variable,levels=tissue_order)
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(tissue_order))
ggplot(to_plot_box, aes(x = cluster, y = value, fill = variable)) +  
  geom_boxplot(outliers = F) + 
  scale_fill_manual(values = color) +
  theme_minimal() +   
  ylab("DNA methylation")+
  ggtitle("Young")

to_plot_box <- tissue_summary
to_plot_box <- to_plot_box[!grepl("chrY", to_plot_box$label), ]
to_plot_box <- merge(to_plot_box,annotation,by="label")
rownames(to_plot_box) <- to_plot_box$label
to_plot_box <- to_plot_box[,-1]
to_plot_box <- reshape2::melt(to_plot_box)
ggplot(to_plot_box, aes(x = cluster, y = value, fill = cluster)) +  
  geom_boxplot(outliers = F) + 
  theme_minimal() +   
  ylab("DNA methylation")+
  ggtitle("Young")

tissue_order <- c("Kidney","Muscle","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung",
                  "Mammary Gland","Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
to_plot_box$variable <- factor(to_plot_box$variable,levels=tissue_order)
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(tissue_order))
ggplot(to_plot_box, aes(x = cluster, y = value, fill = variable)) +  
  geom_boxplot(outliers = F) + 
  scale_fill_manual(values = color) +
  theme_minimal() +   
  ylab("DNA methylation")+
  ggtitle("Young")
