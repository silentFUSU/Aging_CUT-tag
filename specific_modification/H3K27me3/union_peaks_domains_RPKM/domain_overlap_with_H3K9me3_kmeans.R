rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)
library(dplyr)
library(dbplyr)
library(clusterProfiler)
library(GSVA)
library(enrichplot)
options(scipen = 0) 
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

antibody <- "H3K27me3"
annotation <- read.csv("data/samples/all/H3K27me3/edd_domain_merged/kmeans_annotation_RPKM.csv")
split_names <- strsplit(annotation$X, "[:-]")
H3K27me3_df <- data.frame(
  chr = sapply(split_names, "[", 1),
  start = as.numeric(sapply(split_names, "[", 2)),
  end = as.numeric(sapply(split_names, "[", 3)),
  H3K27me3_cluster = annotation$cluster
)

H3K27me3_df <- as.data.table(H3K27me3_df)
setDT(H3K27me3_df)
setkey(H3K27me3_df,chr,start,end)

H3K9me3 <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation.csv")
split_names <- strsplit(H3K9me3$X, "[:-]")
H3K9me3_df <- data.frame(
  chr = sapply(split_names, "[", 1),
  start = as.numeric(sapply(split_names, "[", 2)),
  end = as.numeric(sapply(split_names, "[", 3)),
  H3K9me3_cluster = H3K9me3$cluster
)
H3K9me3_df <- as.data.table(H3K9me3_df)
setDT(H3K9me3_df)
setkey(H3K9me3_df,chr,start,end)
overlaps <- foverlaps(H3K27me3_df, H3K9me3_df, type = "any", nomatch = 0L)  
overlaps <- as.data.frame(overlaps)
calculate_intersection_length <- function(start1, end1, start2, end2) {
  intersection_start <- max(start1, start2)
  intersection_end <- min(end1, end2)
  if (intersection_start < intersection_end) {
    return(intersection_end - intersection_start)
  } else {
    return(0)
  }
}
overlaps$intersection_length <- mapply(calculate_intersection_length, 
                                       overlaps$start, overlaps$end, 
                                       overlaps$i.start, overlaps$i.end)
result <- overlaps %>%
  group_by(H3K27me3_cluster,H3K9me3_cluster) %>%
  summarise(total_intersection_length = sum(intersection_length))

H3K27me3_df <- as.data.frame(H3K27me3_df)
H3K27me3_df$length <- H3K27me3_df$end - H3K27me3_df$start +1
H3K27me3_length <- H3K27me3_df %>%
  group_by(H3K27me3_cluster) %>%
  summarise(total_length = sum(length))

to_plot <- merge(result,H3K27me3_length,by="H3K27me3_cluster")
to_plot$proportion <- to_plot$total_intersection_length/to_plot$total_length * 100

for(i in c(1:3)){
  total_intersection_length <- H3K27me3_length$total_length[which(H3K27me3_length$H3K27me3_cluster==i)] - sum(to_plot$total_intersection_length[which(to_plot$H3K27me3_cluster==i)])
  total_length <- H3K27me3_length$total_length[which(H3K27me3_length$H3K27me3_cluster==i)]
  proportion <- total_intersection_length/total_length*100
  t_to_plot_other <- data.frame(H3K27me3_cluster=i,
                                H3K9me3_cluster="other",
                                total_intersection_length=total_intersection_length,
                                total_length=total_length,
                                proportion=proportion)
  to_plot <- rbind(to_plot,t_to_plot_other)
}
to_plot$H3K27me3_cluster <- paste0("kmeans",to_plot$H3K27me3_cluster)
to_plot$H3K9me3_cluster <- paste0("kmeans",to_plot$H3K9me3_cluster)
to_plot$H3K9me3_cluster[which(to_plot$H3K9me3_cluster=="kmeansother")] <- "other"
color <- setNames(c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF","grey"),c("kmeans1","kmeans2","kmeans3","kmeans4","other"))
ggplot(to_plot, aes(x = H3K27me3_cluster, y = proportion, fill = H3K9me3_cluster)) +  
  geom_bar(stat = 'identity',color="white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20)) +
  ylab("Proportion")+
  xlab("H3K27me3 kmeans cluster")
