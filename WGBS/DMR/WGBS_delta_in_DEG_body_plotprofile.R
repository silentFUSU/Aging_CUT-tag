rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(deepToolsDownstream)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)  
tissue <- args[1]  
condition <- args[2]
se <- importCount(paste0("result/WGBS/",tissue,"/WGBS_change_in_DEG/matrix/WGBS_1kb_delta_in_DEG_body_",condition,".mat.gz"))
if(condition=="increase"){
  color <- setNames("red",tissue)
}else{
  color <- setNames("blue",tissue)
}
se@metadata$sample_labels <-c(tissue)
p <- plotProfile(se,facet = NULL) +scale_x_continuous(  
  labels = c("-10000", "TSS", "TES", "10000")) +
  scale_color_manual(values =  color )+
  theme(text = element_text(size = 15)) +
  ylab("Delta")+
  geom_hline(yintercept = 0, linetype="dashed", color = "red")+
  ggtitle(paste0(tissue,"\nCG% Delta in gene expression ",condition," region"))

ggsave(paste0("result/WGBS/",tissue,"/WGBS_change_in_DEG/plot/WGBS_1kb_delta_in_DEG_body_",condition,".png"),p,width = 8,height = 6,type="cairo")
