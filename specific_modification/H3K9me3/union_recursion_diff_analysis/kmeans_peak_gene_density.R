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
library(GenomicRanges)
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

genes <- fread(paste0("~/ref_data/for_normal_mapping/TSS/refBed/mm10_refGene.bed"))
setDT(genes)
setkey(genes,V1,V2,V3)

### gene number
summary_enrichment <- data.frame()
for(cluster in sort(unique(regions$cluster))){
  df <- regions[which(regions$cluster==cluster),]
  df$label <- paste0(df$chr,":",df$start,"-",df$end)
  length <- data.frame(label=df$label,length=df$end-df$start+1)
  df <- as.data.table(df)
  setDT(df)
  setkey(df,chr,start,end)
  overlaps <- foverlaps(df, genes, type = "any", nomatch = 0L) 
  overlaps <- as.data.frame(overlaps)
  t_summary_enrichment <- as.data.frame(table(overlaps$label))
  colnames(t_summary_enrichment)[1] <- "label"
  t_summary_enrichment <- merge(t_summary_enrichment,length,by="label",all=T)
  t_summary_enrichment$condition <- cluster
  t_summary_enrichment$Freq[is.na(t_summary_enrichment$Freq)] <- 0
  t_summary_enrichment$enrichment_per_kb <- t_summary_enrichment$Freq / t_summary_enrichment$length * 1000
  summary_enrichment <- rbind(summary_enrichment,t_summary_enrichment[,c("condition","enrichment_per_kb")])
}
to_plot <- summary_enrichment
to_plot$condition <- factor(to_plot$condition, levels=c("kmeans1-chrY","kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
color <- setNames(c("#e64b35","#4dbbd5","#00a087","#3c5488","#f39b7f","grey","#636363"),c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
p <- ggplot(to_plot, aes(x = condition, y = enrichment_per_kb,fill=condition)) +
  geom_boxplot(outliers = F) +
  labs(x = NULL, y = "Gene number per kb", title = "Gene density") +
  scale_fill_manual(values = color)+
  theme_bw()+
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
    y_position = c(0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22),
    tip_length = c(1/100,1/100,1/100,1/100,1/100,1/100)
  )
ggsave("result/Sup_figures/H3K9me3_kmeans_gene_number.pdf",p,height = 4,width = 6)

### gene length
summary_percentage <- data.frame()
genes_body <- GRanges(seqnames = genes$V1,ranges = IRanges(start=genes$V2, end=genes$V3))
for(cluster in sort(unique(regions$cluster))){
  df <- regions[which(regions$cluster==cluster),]
  df$length <- df$end - df$start + 1
  df$overlap <- 0
  for(i in c(1:nrow(df))){
    target_peak <- GRanges(seqnames = df$chr[i], 
                           ranges = IRanges(start = df$start[i], end = df$end[i]))
    overlaps <- findOverlaps(target_peak, genes_body)
    overlap_indices <- subjectHits(overlaps)
    overlapping_peaks <- genes_body[overlap_indices]
    target_peak_rep <- rep(target_peak, length(overlapping_peaks))
    intersection_ranges <- pintersect(target_peak_rep, overlapping_peaks)
    merged_ranges <- reduce(intersection_ranges)
    covered_length <- sum(width(merged_ranges))
    df$overlap[i] <- covered_length
  }
  df$gene_overlap_percentage <- df$overlap/df$length * 100
  summary_percentage <- rbind(summary_percentage,df[,c("cluster","gene_overlap_percentage")])
}
to_plot <- summary_percentage
to_plot$cluster <- factor(to_plot$cluster, levels=c("kmeans1-chrY","kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
p <- ggplot(to_plot, aes(x = cluster, y = gene_overlap_percentage,fill=cluster)) +
  geom_boxplot() +
  scale_fill_manual(values = color) +
  labs(x = NULL, y = "Gene length overlap percentage", title = "Gene length overlap percentage") +
  theme_bw()+
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
    y_position = c(110, 120, 130, 140, 150, 160,170),
    tip_length = c(1/100,1/100,1/100,1/100,1/100,1/100)
  )
ggsave("result/Sup_figures/H3K9me3_kmeans_gene_length.pdf",p,height = 4,width = 6)
