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
library(DSS)
library(patchwork)

tissue <- "lung"
DMR <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta01_minCG5.txt"),header = T)
DMR <- DMR[,c(1:3,9)]
DMR$label <- paste0(DMR$chr,":",DMR$start,"-",DMR$end)
DMR_increase <- as.data.table(DMR[which(DMR$areaStat > 0),])
setDT(DMR_increase) 
setkey(DMR_increase,chr,start,end)

DMR_LMJ <- read.table(paste0("data/samples/WGBS/LMJ_DMR_bed/diff_sig_",tissue,"_DMRs.txt"),header = T)
DMR_LMJ$label_LMJ <-paste0(DMR_LMJ$chr,":",DMR_LMJ$start,"-",DMR_LMJ$end)
DMR_increase_LMJ <- DMR_LMJ[which(DMR_LMJ$meth.diff > 0),]
DMR_increase_LMJ <- as.data.table(DMR_increase_LMJ)
setDT(DMR_increase_LMJ)
setkey(DMR_increase_LMJ,chr,start,end)

DMR_increase_overlap <- as.data.frame(foverlaps(DMR_increase_LMJ, DMR_increase, type = "any", nomatch = 0L))
DMR_increase_to_plot <- DMR_increase_overlap[,c("label","areaStat","qvalue")]
DMR_increase_to_plot$areaStat <- log10(abs(DMR_increase_to_plot$areaStat))
DMR_increase_to_plot$qvalue <- -log10(DMR_increase_to_plot$qvalue)

DMR_increase <- as.data.frame(DMR_increase)
DMR_increase_other <- DMR_increase[which(!DMR_increase$label %in% DMR_increase_to_plot$label),c("label","areaStat")]
DMR_increase_other$areaStat <- log10(DMR_increase_other$areaStat)
DMR_increase_other$qvalue <- 0

DMR_increase_LMJ <- as.data.frame(DMR_increase_LMJ)
DMR_increase_LMJ_other <- DMR_increase_LMJ[which(!DMR_increase_LMJ$label_LMJ %in% DMR_increase_overlap$label_LMJ),c("label_LMJ","qvalue")]
DMR_increase_LMJ_other$qvalue <- -log10(DMR_increase_LMJ_other$qvalue)
DMR_increase_LMJ_other$areaStat <- 0

DMR_increase_LMJ_other <- DMR_increase_LMJ_other[,c("label_LMJ","areaStat","qvalue")]
colnames(DMR_increase_LMJ_other)[1] <- "label"
DMR_increase_to_plot <- rbind(DMR_increase_to_plot,DMR_increase_other)
DMR_increase_to_plot <- rbind(DMR_increase_to_plot,DMR_increase_LMJ_other)
ggplot(data = DMR_increase_to_plot, aes(x = abs(areaStat), y = qvalue)) +
  geom_point() +
  labs(
    x = "log10(DSS AreaStat)",
    y = "-log10(methylkit Qvalue)",
    title = "Scatter Plot of AreaStat vs Qvalue"
  ) +
  theme_minimal()
smoothScatter(DMR_increase_to_plot$qvalue ~ DMR_increase_to_plot$areaStat,
              bandwidth = 0.05,ylab = "-log10(methylkit Qvalue)",xlab = "log10(DSS AreaStat)")

DMR_decrease <- as.data.table(DMR[which(DMR$areaStat < 0),])
DMR_decrease$areaStat <- abs(DMR_decrease$areaStat)
setDT(DMR_decrease)
setkey(DMR_decrease,chr,start,end)

DMR_decrease_LMJ <- DMR_LMJ[which(DMR_LMJ$meth.diff < 0),]
DMR_decrease_LMJ <- as.data.table(DMR_decrease_LMJ)
setDT(DMR_decrease_LMJ)
setkey(DMR_decrease_LMJ,chr,start,end)

DMR_decrease_overlap <- as.data.frame(foverlaps(DMR_decrease_LMJ, DMR_decrease, type = "any", nomatch = 0L))
DMR_decrease_to_plot <- DMR_decrease_overlap[,c("label","areaStat","qvalue")]
DMR_decrease_to_plot$areaStat <- log10(abs(DMR_decrease_to_plot$areaStat))
DMR_decrease_to_plot$qvalue <- -log10(DMR_decrease_to_plot$qvalue)

DMR_decrease <- as.data.frame(DMR_decrease)
DMR_decrease_other <- DMR_decrease[which(!DMR_decrease$label %in% DMR_decrease_to_plot$label),c("label","areaStat")]
DMR_decrease_other$areaStat <- log10(DMR_decrease_other$areaStat)
DMR_decrease_other$qvalue <- 0

DMR_decrease_LMJ <- as.data.frame(DMR_decrease_LMJ)
DMR_decrease_LMJ_other <- DMR_decrease_LMJ[which(!DMR_decrease_LMJ$label_LMJ %in% DMR_decrease_overlap$label_LMJ),c("label_LMJ","qvalue")]
DMR_decrease_LMJ_other$qvalue <- -log10(DMR_decrease_LMJ_other$qvalue)
DMR_decrease_LMJ_other$areaStat <- 0

DMR_decrease_LMJ_other <- DMR_decrease_LMJ_other[,c("label_LMJ","areaStat","qvalue")]
colnames(DMR_decrease_LMJ_other)[1] <- "label"
DMR_decrease_to_plot <- rbind(DMR_decrease_to_plot,DMR_decrease_other)
DMR_decrease_to_plot <- rbind(DMR_decrease_to_plot,DMR_decrease_LMJ_other)
ggplot(data = DMR_decrease_to_plot, aes(x = abs(areaStat), y = qvalue)) +
  geom_point() +
  labs(
    x = "log10(DSS AreaStat)",
    y = "-log10(methylkit Qvalue)",
    title = "Scatter Plot of AreaStat vs Qvalue"
  ) +
  theme_minimal()
smoothScatter(DMR_decrease_to_plot$qvalue ~ DMR_decrease_to_plot$areaStat,
              bandwidth = 0.05,ylab = "-log10(methylkit Qvalue)",xlab = "log10(DSS AreaStat)")
