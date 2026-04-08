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

calculate_CpG_content <- function(seq) {
  seq_string <- as.character(seq)  # 转换为字符串
  total_CG <- sum(vcountPattern("CG", seq_string))  # 统计CG出现的次数
  total_bases <- nchar(seq_string)  # 统计总碱基数
  cpg_content <- total_CG / total_bases * 100  # 计算CpG含量
  return(c(total_CG, total_bases, cpg_content))
}
genome <- BSgenome.Mmusculus.UCSC.mm10
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")


df <- read.table("~/ref_data/mm10_10kb_bins.bed")
df$V2 <- df$V2+1
gr <- GRanges(seqnames = df$V1,
              ranges = IRanges(start = df$V2, 
                               end = df$V3))
seqs <- getSeq(genome, gr)
cpg_results <- as.data.frame(t(sapply(seqs, calculate_CpG_content)))
colnames(cpg_results) <- c("CpG_count", "Total_bases", "CpG_content")
rownames(cpg_results) <- df$V4
# cpg_results$cluster <- cluster
summary <- data.frame()
for(tissue in tissues){
  diff <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  diff <- diff[,c("Geneid","Chr","Start","End","Significant")]
  diff <- merge(diff,cpg_results,by.x="Geneid",by.y="row.names")
  diff$tissue <- tissue_label_change(tissue)
  summary <- rbind(summary,diff)
}
saveRDS(summary,"tmp_H3K27me3_CpG_density.rds")

summary <- readRDS("tmp_H3K27me3_CpG_density.rds")
to_plot <- summary[which(summary$tissue=="Lung"),]

ggplot(to_plot, aes(x = Significant, y = CpG_content,fill=Significant)) +
  geom_boxplot(outliers = F) +
  # scale_fill_manual(values = color) +
  labs(x = NULL, y = "CpG_percentage", title = "Lung CpG density") +
  theme_bw()
to_plot <- summary
to_plot$Significant <- factor(to_plot$Significant,levels=c("Up","Stable","Down"))
p<-ggplot(to_plot, aes(x = Significant, y = CpG_content,fill=Significant)) +
  geom_boxplot(outliers = F) +
  # scale_fill_manual(values = color) +
  labs(x = NULL, y = "CpG_percentage", title = "all tissues CpG density") +
  theme_bw()+
  ylim(0,2.5)
p
wilcox.test(to_plot$CpG_content[which(to_plot$Significant=="Up")],to_plot$CpG_content[which(to_plot$Significant=="Down")])
wilcox.test(to_plot$CpG_content[which(to_plot$Significant=="Stable")],to_plot$CpG_content[which(to_plot$Significant=="Down")])

ggsave("result/Sup_figures/all_tissues_H3K27me3_bin_CpG_density.pdf",p,width = 6,height = 8)

#### overlap with H3K27me3 peaks
summary_in_peaks <- data.frame()
for(tissue in tissues){
  t_summary <- summary[which(summary$tissue==tissue_label_change(tissue)),]
  peaks <- read.table(paste0("data/samples/",tissue,"/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed"))
  peaks <- as.data.table(peaks)
  setDT(peaks)
  setkey(peaks,V1,V2,V3)
  bin <- as.data.table(t_summary[,c("Geneid","Chr","Start","End")])
  setDT(bin)
  setkey(bin,Chr,Start,End)
  overlaps <- foverlaps(peaks, bin, type = "any", nomatch = 0L)
  t_summary <- t_summary[which(t_summary$Geneid %in% overlaps$Geneid),]
  summary_in_peaks <- rbind(summary_in_peaks,t_summary)
  }
to_plot <- summary_in_peaks
to_plot$Significant <- factor(to_plot$Significant,levels=c("Up","Stable","Down"))
p<-ggplot(to_plot, aes(x = Significant, y = CpG_content,fill=Significant)) +
  geom_boxplot(outliers = F) +
  # scale_fill_manual(values = color) +
  labs(x = NULL, y = "CpG_percentage", title = "all tissues CpG density") +
  theme_bw()+
  ylim(0,2.5)
p
wilcox.test(to_plot$CpG_content[which(to_plot$Significant=="Up")],to_plot$CpG_content[which(to_plot$Significant=="Down")])
wilcox.test(to_plot$CpG_content[which(to_plot$Significant=="Stable")],to_plot$CpG_content[which(to_plot$Significant=="Down")])

ggsave("result/Sup_figures/all_tissues_H3K27me3_bin_in_young_peaks_CpG_density.pdf",p,width = 6,height = 8)


