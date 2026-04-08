rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(deepToolsDownstream)
args <- commandArgs(trailingOnly = TRUE)  
if (length(args) < 1) {  
  stop("No tissue argument provided")  
}  
file <- args[1]  

se <- importCount(file)
names(se@assays@data@listData) <-c("young1","young2","old1","old2")
se@metadata$sample_labels <-c("young1","young2","old1","old2")
cols  <- setNames(c("red","red","blue","blue"),c("old1","old2","young1","young2"))
plotProfile(se) +scale_x_continuous(  
  labels = c("-10000", "start", "end", "10000")  # 设置标签  
) +scale_color_manual(values = cols) +
  theme(text = element_text(size = 20))
