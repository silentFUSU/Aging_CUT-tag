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
p_list <- list()
conditions <- c("MEF_Bmi1", "MEF_Cbx2", "MEF_Cbx7")
for(condition in conditions){
  H3K27me3 <- read.csv(paste0("data/samples/MEF_OE/H3K27me3/",condition,"/H3K27me3_",condition,"_10kb_bins_diff_after_remove_batch_effect.csv"))
  H2AK119ub1 <- read.csv(paste0("data/samples/MEF_OE/H2AK119ub1/",condition,"/H2AK119ub1_",condition,"_10kb_bins_diff_after_remove_batch_effect.csv"))
  H3K27me3 <- H3K27me3[,c("Geneid","LogFC.oe.vec","Significant")]
  H2AK119ub1 <- H2AK119ub1[,c("Geneid","LogFC.oe.vec","Significant")]
  
  to_plot <- merge(H3K27me3,H2AK119ub1,by="Geneid")
  to_plot$Significant.x[is.na(to_plot$Significant.x)] <- "Stable"
  to_plot$Significant.y[is.na(to_plot$Significant.y)] <- "Stable"
  to_plot$condition <- "Stable"
  to_plot$condition[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Up")] <- "Up"
  to_plot$condition[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Down")] <- "Down"
  to_plot$condition[which(to_plot$Significant.x=="Stable" & to_plot$Significant.y=="Stable")] <- "Stable"
  to_plot$condition[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Down")] <- "Inconsistent"
  to_plot$condition[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Up")] <- "Inconsistent"
  x_range <- range(to_plot$LogFC.oe.vec.x, na.rm = TRUE)  
  y_range <- range(to_plot$LogFC.oe.vec.y, na.rm = TRUE)  
  x_pos_right <- x_range[2] * 0.9    
  x_pos_left <- x_range[1] * 0.9   
  y_pos_top <- y_range[2] * 0.9    
  y_pos_bottom <- y_range[1] * 0.9 
  color <- setNames(c("#e64b35","#3c5488","gray","#00a087"),c("Up","Down","Stable","Inconsistent"))
  p_list[[condition]] <- ggplot(to_plot[which(to_plot$condition=="Stable"),], aes(x = `LogFC.oe.vec.x`, y = `LogFC.oe.vec.y`,color=condition)) +
    geom_point(alpha=0.1) +
    geom_point(data = to_plot[which(to_plot$condition!="Stable"),], aes(x = `LogFC.oe.vec.x`, y = `LogFC.oe.vec.y`,color=condition)) +
    scale_color_manual(values=color)+
    labs(title = paste0(condition),
         x = "H3K27me3",
         y = "H2AK119ub1") +
    theme_bw()+
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Up"),])),  
             x = x_pos_right, y = y_pos_top, colour = "#00b8a9", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Down"),])),  
             x = x_pos_left, y = y_pos_bottom, colour = "#ff9a00", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Up"),])),  
             x = x_pos_left, y = y_pos_top, colour = "#f6416c", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Down"),])),  
             x = x_pos_right, y = y_pos_bottom, colour = "#48466d", size = 5) 
}

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_p <- plot_a_list(p_list,no_of_rows=1,no_of_cols=3)
ggsave("result/MEF_OE/H3K27me3_bin_corr_H2AK119ub1.png",combined_p,width=18,height=5)