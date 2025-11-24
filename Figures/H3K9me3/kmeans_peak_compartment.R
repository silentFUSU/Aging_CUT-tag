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

annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
# annotation <- annotation %>%
#   mutate(cluster = ifelse(grepl("chrY", label) & cluster == "kmeans1", "kmeans1-chrY", cluster))
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

to_plot <- tissue_summary
to_plot <- merge(to_plot,annotation,by="label")
write.csv(to_plot,"data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_mean_compartmentB_percentage_detail.csv")
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
# annotation_color <- list(cluster=setNames(c("#f6416c", "#f8f3d4", "#ffde7d", "#00b8a9","grey"),c(paste0("kmeans",1:4),"Stable")))
pheatmap::pheatmap(to_plot,color = color_palette,cluster_cols = F,cluster_rows = F,show_rownames = F,annotation_row = annotation)

annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
# annotation <- annotation %>%
#   mutate(cluster = ifelse(grepl("chrY", label) & cluster == "kmeans1", "kmeans1-chrY", cluster))
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


tissue_label <- c("Kidney","Muscle","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung",
                  "Mammary Gland","Pancreas","Bone Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(tissue_label))

to_plot_box <- tissue_summary
to_plot_box <- merge(to_plot_box,annotation,by="label")
rownames(to_plot_box) <- to_plot_box$label
to_plot_box <- to_plot_box[,-1]
to_plot_box <- reshape2::melt(to_plot_box)
to_plot_box$cluster <- factor(to_plot_box$cluster,levels=c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
p <-ggplot(to_plot_box, aes(x = cluster, y = value, fill = cluster)) +  
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

to_plot_bar <- tissue_summary
to_plot_bar <- merge(to_plot_bar,annotation,by="label")
rownames(to_plot_bar) <- to_plot_bar$label
# to_plot_bar <- to_plot_bar[,-1]
to_plot_bar <- reshape2::melt(to_plot_bar)
to_plot_bar_tissue <- to_plot_bar %>%
  group_by(cluster, variable) %>%
  summarize(mean_value = mean(value, na.rm = TRUE))
to_plot_bar_tissue_cluster <- to_plot_bar_tissue %>%
  group_by(cluster) %>%
  summarize(mean_value = mean(mean_value, na.rm = TRUE))
write.csv(to_plot_bar_tissue_cluster,paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_mean_compartmentB_percentage.csv"))

to_plot_bar_tissue_cluster$cluster <- factor(to_plot_bar_tissue_cluster$cluster,levels=c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
p <- ggplot(to_plot_bar_tissue_cluster, aes(x = cluster, y = mean_value,fill=cluster)) +
  geom_bar(stat = "identity") +
  geom_point(data = to_plot_bar_tissue, aes(x = cluster, y = mean_value, color = variable), size = 3, position = position_jitter(width = 0.2, height = 0)) +
  scale_color_manual(values = color)+
  theme_bw() +
  labs(x = "Cluster", y = "percentage",title = "compartment B percentage")+
  ylim(0,100)+
  theme(  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14,angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold")
  ) 
ggsave(paste0("result/Sup_figures/H3K9me3_kmeans_compartment.pdf"),p,height = 4,width = 8)
