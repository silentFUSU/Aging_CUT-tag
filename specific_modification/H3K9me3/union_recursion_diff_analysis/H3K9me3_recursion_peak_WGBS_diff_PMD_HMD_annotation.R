rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
library(ggsignif)
library(data.table)
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

calculate_overlap <- function(start1, end1, start2, end2) {
  max_overlap_start <- max(start1, start2)
  min_overlap_end <- min(end1, end2)
  overlap_length <- max(0, min_overlap_end - max_overlap_start)
  return(overlap_length)
}

PMD_HMD_region <- read.table("data/public_data/PMD_coordinates_mm10.bed")
PMD_HMD_region$V2 <- PMD_HMD_region$V2 + 1
PMD_HMD_region$label <- paste0("bin",c(1:nrow(PMD_HMD_region)))
PMD_HMD_region <- PMD_HMD_region[,c("V1","V2","V3","V5","label")]
PMD_HMD_region$V5[is.na(PMD_HMD_region$V5)] <- "other"
PMD_HMD_region <- as.data.table(PMD_HMD_region)
setDT(PMD_HMD_region)
setkey(PMD_HMD_region,V1,V2,V3)
### H3K9me3 peaks annotation
H3K9me3_peaks <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
split_chr <- strsplit(as.character(H3K9me3_peaks$label), ":")  
chr_column <- sapply(split_chr, `[[`, 1)  
split_start_end <- strsplit(sapply(split_chr, `[[`, 2), "-")  
start_column <- sapply(split_start_end, `[[`, 1)  
end_column <- sapply(split_start_end, `[[`, 2)  
H3K9me3_peaks <- data.frame(chr = chr_column,start = start_column, end = end_column, cluster=H3K9me3_peaks$cluster)
H3K9me3_peaks$start <- as.numeric(H3K9me3_peaks$start)
H3K9me3_peaks$end <- as.numeric(H3K9me3_peaks$end)
H3K9me3_peaks <- as.data.table(H3K9me3_peaks)
H3K9me3_peaks$H3K9me3_label <- paste0(H3K9me3_peaks$chr,":",H3K9me3_peaks$start,"-",H3K9me3_peaks$end)
setDT(H3K9me3_peaks)
setkey(H3K9me3_peaks,chr,start,end)
overlaps <- foverlaps(PMD_HMD_region, H3K9me3_peaks, type = "any", nomatch = 0L)  
overlaps <- as.data.frame(overlaps)
result <- overlaps %>%
  group_by(H3K9me3_label) %>%
  summarize(
    total = n(),
    PMD_count = sum(V5 == "PMD"),
    PMD_ratio = PMD_count / total *100
  )
result <- merge(result,H3K9me3_peaks[,c("H3K9me3_label","cluster")],by="H3K9me3_label")
to_plot <- result
to_plot <- to_plot[order(to_plot$cluster),]
ggplot(to_plot, aes(x = cluster, y = PMD_ratio, fill = cluster)) +  
  geom_boxplot(outliers = F) + 
  theme_minimal() +   
  ylab("PMD percentage")+
  ggtitle("Young")


### WGBS changed in H3K9me3 peaks
tissue_mean_summary <- read.csv("data/samples/WGBS/all_tissues_delta_in_100kb_bins_PMD_HMD_cross_comparison_means.csv",row.names = 1)
H3K9me3_peaks <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv")
split_chr <- strsplit(as.character(H3K9me3_peaks$X), ":")  
chr_column <- sapply(split_chr, `[[`, 1)  
split_start_end <- strsplit(sapply(split_chr, `[[`, 2), "-")  
start_column <- sapply(split_start_end, `[[`, 1)  
end_column <- sapply(split_start_end, `[[`, 2)  
H3K9me3_peaks <- data.frame(chr = chr_column,start = start_column, end = end_column, cluster=H3K9me3_peaks$cluster)
H3K9me3_peaks$start <- as.numeric(H3K9me3_peaks$start)
H3K9me3_peaks$end <- as.numeric(H3K9me3_peaks$end)
H3K9me3_peaks <- as.data.table(H3K9me3_peaks)
setDT(H3K9me3_peaks)
setkey(H3K9me3_peaks,chr,start,end)

overlaps <- foverlaps(PMD_HMD_region, H3K9me3_peaks, type = "any", nomatch = 0L)  
overlaps <- as.data.frame(overlaps)
result <- overlaps %>%
  rowwise() %>%
  mutate(overlap_length = calculate_overlap(start, end, V2, V3)) %>%
  ungroup() %>%
  arrange(label, desc(overlap_length)) %>%
  group_by(label) %>%
  slice(1) %>% 
  ungroup()

result$V1 <- factor(result$V1,levels=paste0("chr",c(1:19,"X","Y")))
result$label <- factor(result$label,levels=PMD_HMD_region$label)
result <- result[order(result$cluster, result$V1, result$label), ]

to_plot <- tissue_mean_summary[which(tissue_mean_summary$label %in% result$label),]
to_plot$label <- factor(to_plot$label,levels=result$label)
rownames(to_plot) <- to_plot$label
to_plot <- to_plot[order(to_plot$label),]
to_plot <- to_plot[,-1]

annotation <- as.data.frame(result[,c("label","cluster","V5")])
rownames(annotation) <- annotation$label
annotation <- annotation[,-1]
colnames(annotation) <- c("cluster","condition")
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
annotation_color <- list(cluster=setNames(c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),c(1:4)))
breaks <- c(seq(-0.5, -0.11, length.out = 40), seq(-0.1, 0.1, length.out = 20), seq(0.11, 0.5, length.out = 40))
tissues_order <- c("Kidney","Muscle","Skin","Bladder","Stomach","Heart","Hippocampus","Uterus","Liver","Aorta","Testis","Cortex","Tongue","Cerebellum","BAT","Lung","Mammary.Gland",
                   "Pancreas","Bone.Marrow","iWAT","Cecum","Colon","Jejunum","Spleen","Thymus","Ileum","Ovary")
to_plot <- to_plot[,tissues_order]
annotation$cluster <- factor(annotation$cluster,c(11:4))
pheatmap::pheatmap(to_plot,cluster_rows = F,show_rownames = F,breaks = breaks, annotation_row = annotation,annotation_colors = annotation_color,cluster_cols = F, color = color_palette, clustering_distance_cols="manhattan")


