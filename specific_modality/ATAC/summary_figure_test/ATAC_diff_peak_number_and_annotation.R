rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicRanges)
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peak_num_summary <- list(up=data.frame(),down=data.frame())
peak_annotation_summary <- list(up=data.frame(),down=data.frame())
conditions <- c("up","down")
for(tissue in tissues){
  for(condition in conditions){
    # df <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set/",tissue,".bed"))
    # df <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/bed/ATAC_macs_young_old_narrowpeak_01.bed"))
    df <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set_DAR/",tissue,"_DARs.txt"))
    peaks <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set/",tissue,".bed"))
    rownames(peaks) <- paste0(peaks$V1,":",peaks$V2,"-",peaks$V3)
    df <- merge(df,peaks,by="row.names")
    if(condition == "up"){
      df <- df[which(df$logFC >0 & df$FDR <0.1),]
    }else{
      df <- df[which(df$logFC <0 & df$FDR <0.1),]
    }
    if(nrow(df) >= 20){
      gr <- GRanges(seqnames = df$V1,
                    ranges = IRanges(start = df$V2, end = df$V3))
      peakAnno <- annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db")
      t_peak_num_summary <- data.frame(tissue=tissue_label_change(tissue),num=peakAnno@peakNum)
      peak_num_summary[[condition]] <- rbind(peak_num_summary[[condition]],t_peak_num_summary)
      t_peak_annotation_summary <- peakAnno@annoStat
      colnames(t_peak_annotation_summary)[2] <- tissue_label_change(tissue)
      if(nrow(peak_annotation_summary[[condition]])==0){
        peak_annotation_summary[[condition]] <- t_peak_annotation_summary
      }else{
        peak_annotation_summary[[condition]] <- merge(peak_annotation_summary[[condition]],t_peak_annotation_summary,by="Feature",all=T)
      }
    }
  } 
}

for(condition in conditions){
  peak_annotation_summary[[condition]][is.na(peak_annotation_summary[[condition]])] <- 0
  peak_annotation_summary[[condition]] <- reshape2::melt(peak_annotation_summary[[condition]])
  peak_annotation_summary[[condition]]$condition <- condition
  peak_annotation_summary[[condition]] <- merge(peak_annotation_summary[[condition]], peak_num_summary[[condition]],by.x="variable",by.y="tissue")
  peak_annotation_summary[[condition]]$count <- peak_annotation_summary[[condition]]$value * peak_annotation_summary[[condition]]$num /100
  if(condition=="down"){
    peak_annotation_summary[[condition]]$count <- -peak_annotation_summary[[condition]]$count    
  }
}
num_summary <- merge(peak_num_summary[["up"]],peak_num_summary[["down"]],by="tissue",all=T)
num_summary[is.na(num_summary)] <- 0
num_summary$sum <- rowSums(num_summary[,-1])
num_summary <- num_summary[order(num_summary$sum,decreasing = T),]
to_plot <- rbind(peak_annotation_summary[["up"]],peak_annotation_summary[["down"]])
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,sort(unique(to_plot$Feature)))
to_plot$variable <- factor(to_plot$variable,levels=num_summary$tissue)
ggplot(to_plot, aes(x = variable, y = count, fill = Feature)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme_minimal() +
  scale_fill_manual(values=color) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Feature Count", title = "Peak Annotation Summary by Tissue")
