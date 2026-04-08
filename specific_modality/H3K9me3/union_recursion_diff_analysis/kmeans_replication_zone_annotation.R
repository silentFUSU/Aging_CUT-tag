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

annotation <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
split_names <- strsplit(annotation$label, "[:-]")
annotation_df <- data.frame(
  chr = sapply(split_names, "[", 1),
  start = sapply(split_names, "[", 2),
  end = sapply(split_names, "[", 3),
  cluster = annotation$cluster
)
random_regions1 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed")
colnames(random_regions1) <- c("chr","start","end")
random_regions1$cluster <- "random whole genome"
annotation_df <- rbind(annotation_df,random_regions1)

random_regions2 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed")
colnames(random_regions2) <- c("chr","start","end")
random_regions2$cluster <- "random out of peak"
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
zone_bed <- read.csv("data/public_data/mESC_F121-9_musAllele_IZ.csv",sep = "\t")
colnames(zone_bed) <- c("V1","V2","V3","cluster")
zone_bed<- as.data.table(zone_bed)
setDT(zone_bed)
setkey(zone_bed,V1,V2,V3)
overlaps <- foverlaps(zone_bed, annotation_df, type = "any", nomatch = 0L)  
overlaps <- as.data.frame(overlaps)
result <- overlaps %>%
  group_by(cluster) %>%
  count(i.cluster) %>%
  mutate(percentage = n / sum(n) * 100)

to_plot <- as.data.frame(result)
color <- setNames(c("#e84545","#ff9999","#112d4e","#3f72af"),c("early","earlymid","late","latemid"))
ggplot(to_plot, aes(x = cluster, y = percentage, fill = i.cluster)) +  
  geom_bar(stat = 'identity',color="black") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")



