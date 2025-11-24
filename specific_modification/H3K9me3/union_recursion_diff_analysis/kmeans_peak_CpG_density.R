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

### CpG island
CpG <- read.table("~/ref_data/for_normal_mapping/mm10/cpgi.mm10.bed.txt")
CpG <- CpG[which(CpG$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
CpG <- as.data.table(CpG)
setDT(CpG)
setkey(CpG,V1,V2,V3)
summary_percentage <- data.frame()
for(cluster in sort(unique(regions$cluster))){
  df <- regions[which(regions$cluster==cluster),]
  df <- df[which(df$chr %in% paste0("chr",c(1:19,"X"))),]
  df$label <- paste0(df$chr,":",df$start,"-",df$end)
  length <- data.frame(label=df$label,length=df$end-df$start+1)
  df <- as.data.table(df)
  setDT(df)
  setkey(df,chr,start,end)
  overlaps <- foverlaps(df, CpG, type = "any", nomatch = 0L) 
  overlaps <- as.data.frame(overlaps)
  
  overlaps$CpG_length <- overlaps$V3 - overlaps$V2 +1

  t_summary_percentage <- overlaps %>%
    group_by(label) %>%
    summarize(total_CpG_length = sum(CpG_length, na.rm = TRUE))
  
  t_summary_percentage <- merge(t_summary_percentage,length,by="label",all=T)
  t_summary_percentage$condition <- cluster
  t_summary_percentage$total_CpG_length[is.na(t_summary_percentage$total_CpG_length)] <- 0 
  t_summary_percentage$CpG_length_percentage <- t_summary_percentage$total_CpG_length/t_summary_percentage$length*100
  summary_percentage <- rbind(summary_percentage,t_summary_percentage[,c("condition","CpG_length_percentage")])
}
to_plot <- summary_percentage
to_plot$condition <- factor(to_plot$condition,levels=c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
ggplot(to_plot, aes(x = condition, y = CpG_length_percentage,fill=condition)) +
  geom_boxplot(outliers = F) +
  labs(x = NULL, y = "CpG_island_length_percentage", title = "CpG island density") +
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

#### CG density
calculate_CpG_content <- function(seq) {
  seq_string <- as.character(seq)  # 转换为字符串
  total_CG <- sum(vcountPattern("CG", seq_string))  # 统计CG出现的次数
  total_bases <- nchar(seq_string)  # 统计总碱基数
  cpg_content <- total_CG / total_bases * 100  # 计算CpG含量
  return(c(total_CG, total_bases, cpg_content))
}
genome <- BSgenome.Mmusculus.UCSC.mm10
# regions$cluster[which(regions$chr == "chrY" & regions$cluster == "kmeans1")] <- "kmeans1-chrY"
summary <- data.frame()
for(cluster in sort(unique(regions$cluster))){
  df <- regions[which(regions$cluster==cluster),]
  gr <- GRanges(seqnames = df$chr,
                ranges = IRanges(start = df$start, 
                                 end = df$end))
  seqs <- getSeq(genome, gr)
  cpg_results <- as.data.frame(t(sapply(seqs, calculate_CpG_content)))
  colnames(cpg_results) <- c("CpG_count", "Total_bases", "CpG_content")
  cpg_results$cluster <- cluster
  summary <- rbind(summary,cpg_results)    
}
to_plot <- summary
to_plot$cluster <- factor(to_plot$cluster,levels=c("kmeans1-chrY","kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))
color <- setNames(c("#e64b35","#4dbbd5","#00a087","#3c5488","#f39b7f","grey","#636363"),c("kmeans1","kmeans2","kmeans3","kmeans4","Stable","random out of peak","random whole genome"))

p <- ggplot(to_plot, aes(x = cluster, y = CpG_content,fill=cluster)) +
  geom_boxplot(outliers = F) +
  scale_fill_manual(values = color) +
  labs(x = NULL, y = "CpG_percentage", title = "CpG density") +
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
    y_position = c(1, 1.2, 1.4,1.6,1.8, 2, 2.2),
    tip_length = c(1/100,1/100,1/100,1/100,1/100,1/100)
  )
ggsave(paste0("result/Sup_figures/H3K9me3_kmeans_CpG_density.pdf"),p,height = 4,width = 6)

