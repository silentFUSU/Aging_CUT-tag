rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(factoextra)
library(cluster)
library(umap)
options(bitmapType = "cairo")  
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
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))
tissue_summary <- data.frame()
annotation_col <- data.frame()
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
  tab = read.delim(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak_summits_spm3_all_tissues_merge.counts"),skip=1)
  if(tissue %in% c("mammarygland","uterus","ovary")){
    tab <- tab[which(tab$Chr %in% c(paste0("chr",c(1:19,"X")))),]
  }
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  counts <- counts[,search_table$sample_name]
  
  CPM <- as.data.frame(edgeR::cpm(counts))
  young_cpm <- CPM[,search_table$sample_name[which(search_table$age=="3m")]]
  old_cpm <- CPM[,search_table$sample_name[which(search_table$age=="24m")]]
  young_cpm$rowmeans <- rowMeans(young_cpm) 
  old_cpm$rowmeans <- rowMeans(old_cpm)
  t_tissue_summary <- merge(young_cpm[,"rowmeans",drop=F],old_cpm[,"rowmeans",drop=F],by="row.names",all=T)
  colnames(t_tissue_summary) <- c("label",paste0(tissue_label_change(tissue),"_young"),paste0(tissue_label_change(tissue),"_old"))
  if(nrow(tissue_summary)==0){
    tissue_summary <- t_tissue_summary
  }else{
    tissue_summary <- merge(tissue_summary,t_tissue_summary,by="label")
  }
  t_annotation_col <- data.frame(sample=c(paste0(tissue_label_change(tissue),"_young"), paste0(tissue_label_change(tissue),"_old")),
                                 tissue=c(tissue_label_change(tissue),tissue_label_change(tissue)))
  annotation_col <- rbind(annotation_col,t_annotation_col)
  }

# to_plot <- tissue_summary[c(1:1000),]
to_plot <- tissue_summary
column_names <- c(
  "label", "Aorta_young", "Aorta_old", "BAT_young", "BAT_old",
  "Bladder_young", "Bladder_old", "Bone Marrow_young", "Bone Marrow_old", 
  "Cecum_young", "Cecum_old","Cerebellum_young", "Cerebellum_old", 
  "Colon_young", "Colon_old","Cortex_young","Cortex_old", 
  "Heart_young", "Heart_old", "Hippocampus_young", "Hippocampus_old", 
  "Ileum_young", "Ileum_old", "iWAT_young", "iWAT_old",
  "Jejunum_young", "Jejunum_old", "Kidney_young", "Kidney_old", 
  "Liver_young", "Liver_old", "Lung_young", "Lung_old", 
  "Mammary Gland_young", "Mammary Gland_old","Muscle_young", "Muscle_old", 
  "Ovary_young", "Ovary_old", "Pancreas_young","Pancreas_old", 
  "Skin_young", "Skin_old", "Spleen_young", "Spleen_old",
  "Stomach_young", "Stomach_old", "Testis_young", "Testis_old", 
  "Thymus_young","Thymus_old", "Tongue_young", "Tongue_old", 
  "Uterus_young", "Uterus_old"
)
to_plot <- to_plot[,column_names]
rownames(to_plot) <- to_plot$label
to_plot <- to_plot[,-1]
k <- 28
kmeans_df <- to_plot
set.seed(1)
kmeans_result <- kmeans(kmeans_df, centers=k)
kmeans_result <- as.data.frame(kmeans_result$cluster)
colnames(kmeans_result) <- "cluster"
to_plot <- merge(to_plot,kmeans_result,by="row.names") 
rownames(to_plot) <- to_plot$Row.names
to_plot <- to_plot[,-1]
to_plot <- to_plot[order(to_plot$cluster),]

kmeans_average_value <- data.frame()
for(i in c(1:k)){
  t_kmeans_average_value <- to_plot[which(to_plot$cluster==i),-which(colnames(to_plot) %in% c("cluster"))]
  t_kmeans_average_value <- data.frame(kmeans=i,mean=mean(unlist(t_kmeans_average_value), na.rm = T))
  kmeans_average_value <- rbind(kmeans_average_value,t_kmeans_average_value)
}
kmeans_average_value <- kmeans_average_value[order(kmeans_average_value$mean,decreasing = T),]

high_kmeans <- kmeans_average_value$kmeans[which(kmeans_average_value$mean>6)]
high_kmeans <- c(20, 19,  4,  1, 18, 12, 22, 24,  9, 15, 17,11,25,16,3,10,26,21,8,2,23,27,14)
kmeans_order <- data.frame()
for(i in c(1:k)){
  df <- to_plot[which(to_plot$cluster==i),-which(colnames(to_plot)=="cluster")]
  max_index <- which(df == max(df), arr.ind = TRUE)
  max_col_name <- names(df)[max_index[2]]
  t_kmeans_order <- data.frame(cluster=i,colnames=max_col_name)
  kmeans_order <- rbind(kmeans_order,t_kmeans_order)
}
kmeans_order <- kmeans_order[which(!kmeans_order$cluster %in% high_kmeans),]
kmeans_order <- kmeans_order %>%
  separate(colnames, into = c("tissue", "age_group"), sep = "_")
kmeans_order$tissue <- factor(kmeans_order$tissue,levels = sapply(tissues,tissue_label_change))
kmeans_order <- kmeans_order[order(kmeans_order$tissue),]
kmeans_order <- c(high_kmeans,2,7,6,13,5)
to_plot$cluster <- factor(to_plot$cluster,levels = c(high_kmeans,7,6,13,5))
to_plot <- to_plot[order(to_plot$cluster),]
annotation_row <- to_plot[,c("cluster"),drop=F]
to_plot <- to_plot[,-which(colnames(to_plot)=="cluster")]

annotation_row$cluster <- as.character(annotation_row$cluster)
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- c(seq(0, 3.5, length.out = 40), seq(3.6, 6.5, length.out = 20), seq(6.6, 10, length.out = 40))

color_row <- read.table("data/samples/30_distinct_color.txt")
color_row <- setNames(color_row$V1,as.character(c(1:k)))
color_col <- read.table("data/samples/30_distinct_color.txt")
color_col <- setNames(color_col$V1[1:27],sort(unique(annotation_col$tissue)))
annotation_color <- list(cluster=color_row,tissue=color_col)

# rownames(annotation_col) <- annotation_col$sample
# annotation_col <- annotation_col[,-1,drop=F]

p <- pheatmap::pheatmap(to_plot,show_rownames = F,breaks = breaks, 
                        annotation_row = annotation_row,annotation_col = annotation_col,
                        color = color_palette,cluster_rows = F,cluster_cols = F,
                        annotation_colors = annotation_color,filename = "result/all/ATAC/all_tissues_macs_young_old_narrowpeak_summits_spm3_all_tissues_merge.png",width = 8,height =12)


# p <- pheatmap::pheatmap(to_plot,show_rownames = F,breaks = breaks, 
#                         annotation_row = annotation_row,annotation_col = annotation_col,
#                         color = color_palette,cluster_rows = F,cluster_cols = F,
#                         annotation_colors = annotation_color,filename = "result/all/ATAC/all_tissues_macs_young_old_narrowpeak_summits_spm3_all_tissues_merge.png",width = 40,height =80)
# 
