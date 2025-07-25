rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(universalmotif)

motif <- readRDS("data/public_data/cisBP_mouse_pfms_2021.rds")
write_homer(motif,"~/ref_data/for_normal_mapping/mm10/homer_cisbp_05.motif",overwrite = T,threshold = 0.5)
write_meme(motif,"~/ref_data/for_normal_mapping/mm10/meme_cisbp_05.meme")
