rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(stringr)buion_order)
color <- setNames(c("red","grey","blue"),c("Up","Stable","Down"))
p <- ggplot(summary,mapping = aes(x=`LogFC.oe.vec`,y=condition,fill = Significant))+
  geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+
  scale_fill_manual(values=color)+
  theme(text = element_text(size = 13))+ 
  ggtitle(paste0(gene," RNA")) 
p
ggsave("result/figures/MEF_OE_Cdkn2a_RNA.pdf",p,width = 6,height = 8)
