rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(rtracklayer)
library(gridExtra)
library(grid)  
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggsignif)
options(bitmapType="cairo")  

regions <- read.csv("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv",row.names = 1)
split_chr <- strsplit(as.character(regions$label), ":")  
chr_column <- sapply(split_chr, `[[`, 1)  
split_start_end <- strsplit(sapply(split_chr, `[[`, 2), "-")  
start_column <- sapply(split_start_end, `[[`, 1)  
end_column <- sapply(split_start_end, `[[`, 2)  
regions <- data.frame(chr = chr_column,start = start_column, end = end_column, cluster=regions$cluster)


random_regions1 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY.bed")
random_regions1$cluster <- "random whole genome"
colnames(random_regions1)[1:3] <- c("chr","start","end")
regions <- rbind(regions,random_regions1)

random_regions2 <- read.table("data/samples/all/H3K9me3/bed/H3K9me3_random_200kb_region_rmchrY_out_recursion_peaks.bed")
random_regions2$cluster <- "random out of peak"
colnames(random_regions2)[1:3] <- c("chr","start","end")
regions <- rbind(regions,random_regions2)

regions$start <- as.numeric(regions$start)
regions$end <- as.numeric(regions$end)
regions <- regions[which(regions$chr %in% paste0("chr",c(1:19,"X","Y"))),]
# regions$cluster[which(regions$chr == "chrY" & regions$cluster == "kmeans1")] <- "kmeans1-chrY"
# write.table(regions,"data/samples/all/H3K9me3/recursion_peaks_diff_table/bed/all_regions_add_random.bed",append = F,quote = F,sep = "\t",row.names = F,col.names = F)
phastCons <- fread(paste0("~/ref_data/for_normal_mapping/mm10/phastCons/mm10.60way.phastCons.10kb.bedgraph"))
setDT(phastCons)
setkey(phastCons,V1,V2,V3)
summary <- data.frame()
for(cluster in sort(unique(regions$cluster))){
  df <- as.data.table(regions[which(regions$cluster==cluster),])
  df$label <- paste0(df$chr,":",df$start,"-",df$end)
  setDT(df)
  setkey(df,chr,start,end)
  overlaps <- as.data.frame(foverlaps(df, phastCons, type = "any", nomatch = 0L))
  result <- overlaps %>%
    group_by(label) %>%
    summarise(avg_V4 = mean(V4))
  result$cluster <- cluster
  summary <- rbind(summary,result)
}

to_plot <- as.data.frame(summary)
write.csv(to_plot,"data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY/kmeans_median_phastCons_score_detail.csv")
to_plot$cluster <- factor(to_plot$cluster,levels=c("kmeans1-chrY","kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
ggplot(to_plot, aes(x = cluster, y = avg_V4,fill=cluster)) +
  geom_boxplot(outliers = F) +
  labs(x = NULL, y = "phastCons score", title = "phastCons score") +
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 12,face = "bold", color = "black", angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 12,face = "bold", color = "black"),  
    axis.title.x = element_text(size = 14,face = "bold", color = "black"), 
    axis.title.y = element_text(size = 14,face = "bold", color = "black"), 
    legend.text = element_text(size = 12),  
    panel.background = element_blank(),  
    panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "grey"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "lightgrey"),
    plot.title = element_text(size = 16, face = "bold")
  )+
  geom_signif(
    comparisons = list(c("kmeans1", "kmeans2"), c("kmeans1", "kmeans3"),c("kmeans1", "kmeans4"),c("kmeans1", "Stable"),c("kmeans1", "random out of peak"),c("kmeans1", "random whole genome")),
    textsize = 4,test = "t.test",
    map_signif_level = TRUE,
    y_position = c(1, 1.2, 1.4,1.6,1.8, 2, 2.2),
    tip_length = c(1/100,1/100,1/100,1/100,1/100,1/100)
  )
result <- to_plot %>%
  group_by(cluster) %>%
  summarise(median_score = median(avg_V4, na.rm = TRUE))
write.csv(result,paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/feature_table_combined_chrY//kmeans_median_phastCons_score.csv"))

