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
senescence <- read.csv("data/samples/RNA/MEF_mid_age/diff_expression_gene_strict_filter_bar.csv")
senescence <- senescence[,c("X","logFC","Significant")]
# conditions <- c("MEF_Bmi1", "MEF_Cbx2", "MEF_Cbx7")
# conditions <- c()
conditions <- c("MEF_Bmi1", "MEF_Cbx2", "MEF_Cbx7", "MEF_mCbx8", "MEF_mEzh2", "MEF_hEzh2")
p_list <- list()
cor_summary <- data.frame()
color_df <- data.frame(condition=c("MEF_Bmi1", "MEF_Cbx2", "MEF_Cbx7", "MEF_mCbx8","MEF_mEzh2","MEF_hEzh2"),color=c("#FF4D1A","#9B111E","#4A0713","#E0002A","#3c5488","#4dbbd5"))
for(condition in conditions){
  if(condition %in% c("MEF_mEzh2", "MEF_hEzh2")){
    df <- read.csv(paste0("data/samples/RNA/MEF_OE_RNA/MEF_Vector2_",condition,"_diff_expression_gene_strict_filter_bar.csv"))
  }else{
    df <- read.csv(paste0("data/samples/RNA/MEF_OE_RNA/MEF_Vector_",condition,"_diff_expression_gene_strict_filter_bar.csv"))
  }
  df <- df[,c("Geneid","LogFC.oe.vec","Significant")]
  to_plot <- merge(df,senescence,by.x="Geneid", by.y="X")
  to_plot$Significant.x[is.na(to_plot$Significant.x)] <- "Stable"
  to_plot$Significant.y[is.na(to_plot$Significant.y)] <- "Stable"
  to_plot$condition <- "Stable"
  to_plot$condition[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Up")] <- "Up"
  to_plot$condition[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Down")] <- "Down"
  to_plot$condition[which(to_plot$Significant.x=="Stable" & to_plot$Significant.y=="Stable")] <- "Stable"
  to_plot$condition[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Down")] <- "Inconsistent"
  to_plot$condition[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Up")] <- "Inconsistent"
  x_range <- range(to_plot$LogFC.oe.vec, na.rm = TRUE)  
  y_range <- range(to_plot$logFC, na.rm = TRUE)  
  x_pos_right <- x_range[2] * 0.9    
  x_pos_left <- x_range[1] * 0.9   
  y_pos_top <- y_range[2] * 0.9    
  y_pos_bottom <- y_range[1] * 0.9 
  color <- setNames(c("#e64b35","#3c5488","gray","#00a087"),c("Up","Down","Stable","Inconsistent"))
  
  cor <- cor.test(to_plot$LogFC.oe.vec,to_plot$logFC)
  t_cor_summary <- data.frame(condition=condition,cor=cor$estimate,pvalue=cor$p.value)
  cor_summary <- rbind(cor_summary,t_cor_summary)
  # p_list[[condition]] <- ggplot(to_plot[which(to_plot$condition=="Stable"),], aes(x = `LogFC.oe.vec`, y = `logFC`,color=condition)) +
  #   geom_point(alpha=0.1) +
  #   geom_point(data = to_plot[which(to_plot$condition!="Stable"),], aes(x = `LogFC.oe.vec`, y = `logFC`,color=condition)) +
  #   scale_color_manual(values=color)+
  #   labs(title = paste0(condition," vs ","P10_P6"),
  #        x = paste0(condition," OE"),
  #        y = "P10 vs P6") +
  #   theme_bw()+
  #   annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Up"),])),
  #            x = x_pos_right, y = y_pos_top, colour = "#00b8a9", size = 5) +
  #   annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Down"),])),
  #            x = x_pos_left, y = y_pos_bottom, colour = "#ff9a00", size = 5) +
  #   annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Down" & to_plot$Significant.y=="Up"),])),
  #            x = x_pos_left, y = y_pos_top, colour = "#f6416c", size = 5) +
  #   annotate("text", label = paste0(nrow(to_plot[which(to_plot$Significant.x=="Up" & to_plot$Significant.y=="Down"),])),
  #            x = x_pos_right, y = y_pos_bottom, colour = "#48466d", size = 5)+
  #   xlim(-8,8)+
  #   ylim(-8,8)
  # fit0 <- lm(`logFC` ~ 0 + `LogFC.oe.vec`, data = to_plot)
  # b <- coef(fit0)[1]
  
  p_list[[condition]] <- ggplot(to_plot, aes(x = `LogFC.oe.vec`, y = `logFC`)) +
    geom_point(color=color_df$color[which(color_df$condition==condition)]) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, color = "black") +
    labs(title = paste0(condition," vs ","P10_P6"),
         x = paste0(condition," OE"),
         y = "P10 vs P6") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_x_continuous(limits = c(-6, 6), breaks = c(-5, 0, 5)) +
    scale_y_continuous(limits = c(-8, 8))
  
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_p <- plot_a_list(p_list,no_of_rows=2,no_of_cols=2)
ggsave("result/figures/MEF_OE_RNA_corr_with_senescence_strict_filter_bar3.pdf",combined_p,width=10,height=9.5)
