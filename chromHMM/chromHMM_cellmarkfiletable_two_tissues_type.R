rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
cellmarkfiletable <- read.table("data/samples/all/combined_analysis_enhancer/cellmarkfiletable.txt")
tissue_type1 <- c("liver","bonemarrow","heart","skin","spleen","cecum","colon","lung","brain","Hip","aorta","muscle","stomach")
tissue_type1_label <- paste0(rep(tissue_type1,4),c("_young1","_young2","_old1","_old2"))
cellmarkfiletable_type1 <- cellmarkfiletable[which(cellmarkfiletable$V1 %in% tissue_type1_label),]
write.table(cellmarkfiletable_type1, "data/samples/all/combined_analysis_enhancer/cellmarkfiletable_type1.txt",sep = "\t",append = F,quote = F,row.names = F,col.names = F)

tissue_type2 <- c("ovary","mammarygland","tongue","uterus","thymus","jejunum","testis","iWAT","BAT","kidney","CB","bladder","pancreas")
tissue_type2_label <- paste0(rep(tissue_type2,4),c("_young1","_young2","_old1","_old2"))
cellmarkfiletable_type2 <- cellmarkfiletable[which(cellmarkfiletable$V1 %in% tissue_type2_label),]
write.table(cellmarkfiletable_type2, "data/samples/all/combined_analysis_enhancer/cellmarkfiletable_type2.txt",sep = "\t",append = F,quote = F,row.names = F,col.names = F)
