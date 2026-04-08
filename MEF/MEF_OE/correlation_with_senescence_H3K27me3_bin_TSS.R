rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
mm10 <- read.table("~/ref_data/mm10_10kb_bins.bed")
mm10 <- as.data.table(mm10)
setDT(mm10)
setkey(mm10,V1,V2,V3)

senescence <- read.csv("data/samples/MEF/H3K27me3/H3K27me3_gene_TSS_10kb_diff_after_remove_batch_effect.csv")
senescence <- senescence[,c("Geneid","LogFC.old.young","Significant")]

antibody <- "H3K27me3"
conditions <- c("MEF_Bmi1", "MEF_Cbx2", "MEF_Cbx7")
p_list <- list()
p_list2 <- list()
peaks <- read.table("data/samples/MEF/H3K27me3/bed/H3K27me3_young_merge-W5000-G10000-E100.bed")
# peaks <- read.table("data/samples/MEF_OE/H2AK119ub1/common_peaks.bed")
peaks <- as.data.table(peaks)
setDT(peaks)
setkey(peaks,V1,V2,V3)

domains <- read.table("data/samples/MEF/H3K27me3/peaks/edd/edd_peaks_fdr05.bed")
domains <- as.data.table(domains)
setDT(domains)
setkey(domains,V1,V2,V3)
for(condition in conditions){
  df <- read.csv(paste0("data/samples/MEF_OE/H3K27me3/",condition,"/H3K27me3_",condition,"_gene_TSS_10kb_diff_after_remove_batch_effect.csv"))
  # df <- df[,c("Geneid","LogFC.oe.vec","Significant")]
  to_plot <- merge(df,senescence,by="Geneid")
  to_plot$Significant.x[is.na(to_plot$Significant.x)] <- "Stable"
  to_plot$Significant.y[is.na(to_plot$Significant.y)] <- "Stable"
  to_plot$condition <- "Stable"
  to_plot$condition[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Up")] <- "Up"
  to_plot$condition[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Down")] <- "Down"
  to_plot$condition[which(to_plot$Significant.x=="Stable" & to_plot$Significant.y=="Stable")] <- "Stable"
  to_plot$condition[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Down")] <- "Inconsistent"
  to_plot$condition[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Up")] <- "Inconsistent"
  x_range <- range(to_plot$LogFC.oe.vec, na.rm = TRUE)  
  y_range <- range(to_plot$LogFC.old.young, na.rm = TRUE)  
  x_pos_right <- x_range[2] * 0.9    
  x_pos_left <- x_range[1] * 0.9   
  y_pos_top <- y_range[2] * 0.9    
  y_pos_bottom <- y_range[1] * 0.9 
  color <- setNames(c("#e64b35","#3c5488","gray","#00a087"),c("Up","Down","Stable","Inconsistent"))
  p_list[[condition]] <- ggplot(to_plot[which(to_plot$condition=="Stable"),], aes(x = `LogFC.oe.vec`, y = `LogFC.old.young`,color=condition)) +
    geom_point(alpha=0.1) +
    geom_point(data = to_plot[which(to_plot$condition!="Stable"),], aes(x = `LogFC.oe.vec`, y = `LogFC.old.young`,color=condition)) +
    scale_color_manual(values=color)+
    labs(title = paste0(condition," vs ","Senescence"),
         x = paste0(condition," OE"),
         y = "Senescence") +
    theme_bw()+
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Up"),])),  
             x = x_pos_right, y = y_pos_top, colour = "#00b8a9", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Down"),])),  
             x = x_pos_left, y = y_pos_bottom, colour = "#ff9a00", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Up"),])),  
             x = x_pos_left, y = y_pos_top, colour = "#f6416c", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Down"),])),  
             x = x_pos_right, y = y_pos_bottom, colour = "#48466d", size = 5) 
  
  bins <- as.data.table(to_plot[,c(1:4)])
  setDT(bins)
  setkey(bins,Chr,Start,End)
  overlaps <- as.data.frame(foverlaps(bins,domains, type = "any", nomatch = 0L))
  to_plot$overlap <- "other"
  to_plot$overlap[which(to_plot$Geneid %in% overlaps$Geneid)] <- "domain"
  
  overlaps <- as.data.frame(foverlaps(bins,peaks, type = "any", nomatch = 0L))
  to_plot$overlap[which(to_plot$Geneid %in% overlaps$Geneid)] <- "peak"
  
  to_plot$quadrant <- "other"
  to_plot$quadrant[which(to_plot$condition=="Up")] <- "first"
  to_plot$quadrant[which(to_plot$condition=="Inconsistent" & to_plot$LogFC.oe.vec < 0)] <- "second"
  to_plot$quadrant[which(to_plot$condition=="Down")] <- "third"
  to_plot$quadrant[which(to_plot$condition=="Inconsistent" & to_plot$LogFC.oe.vec > 0)] <- "fourth"
  
  to_plot_table <- as.data.frame(table(to_plot$overlap, to_plot$quadrant),
                                 stringsAsFactors = FALSE)
  to_plot_table$prop <- with(to_plot_table,
                             ave(Freq, Var2, FUN = function(x) x / sum(x)*100))
  
  
  to_plot_table$Var1 <- factor(to_plot_table$Var1,levels=c("domain","peak","other"))
  to_plot_table$Var2 <- factor(to_plot_table$Var2,levels=c("first","second","third","fourth","other"))
  p_list2[[condition]] <- ggplot(to_plot_table, aes(x = Var2, y = prop, fill = Var1)) +
    geom_col(position = "stack", width = 0.8) +
    labs(x = "Var2", y = "Proportion", fill = "Var1") +
    theme_classic()+ggtitle(condition)
}

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_p <- plot_a_list(p_list,no_of_rows=1,no_of_cols=3)
ggsave("result/MEF_OE/H3K27me3_TSS_bin_corr_with_senescence.png",combined_p,width=18,height=5)

combined_p <- plot_a_list(p_list2,no_of_rows=1,no_of_cols=3)
ggsave("result/MEF_OE/H3K27me3_TSS_bin_corr_with_senescence_overlap.png",combined_p,width=18,height=5)

