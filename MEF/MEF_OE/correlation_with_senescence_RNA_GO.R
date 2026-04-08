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
conditions <- c("MEF_Bmi1", "MEF_Cbx2", "MEF_Cbx7","MEF_mEzh2", "MEF_hEzh2", "MEF_mCbx8")
p_list <- list()
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
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
  
  gene_list <- list()
  
  gene_list[["first"]] <- to_plot$Geneid[which(to_plot$condition=="Up")]
  gene_list[["second"]] <- to_plot$Geneid[which(to_plot$condition == "Inconsistent" & to_plot$LogFC.oe.vec < 0)]
  gene_list[["third"]] <- to_plot$Geneid[which(to_plot$condition == "Down")]
  gene_list[["fourth"]] <- to_plot$Geneid[which(to_plot$condition == "Inconsistent" & to_plot$LogFC.oe.vec > 0)]
  quadrants <- c("first","second","third","fourth")
  p_list <- list()
  for(quadrant in quadrants){
    genes <- bitr(gene_list[[quadrant]],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
    genes_GO <- enrichGO( genes$ENTREZID,#GO富集分析
                                  OrgDb = GO_database,
                                  keyType = "ENTREZID",#设定读取的gene ID类型
                                  ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                                  pvalueCutoff = 0.05,#设定p值阈值
                                  qvalueCutoff = 0.05,#设定q值阈值
                                  readable = T)
    if(nrow(genes_GO) >0){
      p_list[[quadrant]] <- barplot(genes_GO,label_format = 50,showCategory = 10)+ggtitle(paste0(condition," ",quadrant))
    }else{
      p_list[[quadrant]] <- ggplot()
    }
  }
  plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
    
    patchwork::wrap_plots(master_list_with_plots, 
                          nrow = no_of_rows, ncol = no_of_cols)
  }
  combined_p <- plot_a_list(p_list,no_of_rows=2,no_of_cols=2)
  ggsave(paste0("result/MEF_OE/",condition,"_RNA_corr_with_senescence_strict_filter_bar_GO.png"),combined_p,width=18,height=10)
  
}

