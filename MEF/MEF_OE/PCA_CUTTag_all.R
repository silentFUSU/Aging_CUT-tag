rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(patchwork)
library(edgeR)
library(MASS) 
library(Seurat)
library(gridExtra)
library(ggrepel)
library(stringr)

antibodys <- c("H3K27me3","H2AK119ub1")

p_list <- list()
for(antibody in antibodys){
  tab <- read.table(paste0("data/samples/MEF_OE/",antibody,"/",antibody,"_all_10kb.counts"),header = T)
  new_names <- sub("^.*\\.bam\\.([A-Za-z]+\\d+\\.\\d+)\\.bam.*$", "\\1", colnames(tab)[7:ncol(tab)])
  new_names <- gsub("\\.", "-", new_names)
  colnames(tab)[7:ncol(tab)] <- new_names
  tab <- tab[!grepl("chrY", tab$Chr), ]
  rownames(tab) <- tab$Geneid
  tab <- tab[,-1]
  counts <- tab[,c(6:ncol(tab))]
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  search_table$age <- c("vec","vec","Bmi1","Bmi1","Cbx2","Cbx2","Cbx7","Cbx7")
  counts <- counts[,search_table$sample_name]
  y= DGEList(counts=counts,group=search_table$age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  
  logCPMs <- cpm(y, log = TRUE)
  pca <- prcomp(t(logCPMs))
  to_plot <- data.frame(pca$x)
  to_plot$sample_name <- rownames(to_plot)
  
  to_plot <- merge(to_plot,search_table,by="sample_name")
  # to_plot$age <- factor(to_plot$age,levels=c("vec","oe"))
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  p_list[[paste0(antibody)]] <- ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
    ggtitle(paste0(antibody))
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_p <- plot_a_list(p_list,no_of_rows=2,no_of_cols=1)
combined_p
