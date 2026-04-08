rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)

peaks_list <- list()
peaks_list[["Vector"]] <- read.table("data/samples/MEF_OE/H2AK119ub1/MEF_Vector/bed/H2AK119ub1_MEF_Vector_merge-W5000-G10000-E100.bed")
peaks_list[["Bmi1"]] <- read.table("data/samples/MEF_OE/H2AK119ub1/MEF_Bmi1/bed/H2AK119ub1_MEF_Bmi1_merge-W5000-G10000-E100.bed")
peaks_list[["Cbx2"]] <- read.table("data/samples/MEF_OE/H2AK119ub1/MEF_Cbx2/bed/H2AK119ub1_MEF_Cbx2_merge-W5000-G10000-E100.bed")
peaks_list[["Cbx7"]] <- read.table("data/samples/MEF_OE/H2AK119ub1/MEF_Cbx7/bed/H2AK119ub1_MEF_Cbx7_merge-W5000-G10000-E100.bed")

peaks_list[["Vector"]] <- as.data.table(peaks_list[["Vector"]])
setDT(peaks_list[["Vector"]])
setkey(peaks_list[["Vector"]],V1,V2,V3)
common_peaks_list <- list()
for(condition in c("Bmi1","Cbx2","Cbx7")){
  peaks_list[[condition]] <- as.data.table(peaks_list[[condition]][,1:3])
  setDT(peaks_list[[condition]])
  setkey(peaks_list[[condition]],V1,V2,V3)
  overlaps <- as.data.frame(foverlaps(peaks_list[["Vector"]],peaks_list[[condition]], type = "any", nomatch = 0L))
  common_peaks_list[[condition]] <- overlaps$V4
}
peaks <- intersect(common_peaks_list[["Bmi1"]],common_peaks_list[["Cbx2"]])
peaks <- intersect(peaks,common_peaks_list[["Cbx7"]])
vector <- read.table("data/samples/MEF_OE/H2AK119ub1/MEF_Vector/bed/H2AK119ub1_MEF_Vector_merge-W5000-G10000-E100.bed")
peaks <- vector[which(vector$V4 %in% peaks),]
write.table(peaks,"data/samples/MEF_OE/H2AK119ub1/common_peaks.bed",append = F,quote = F,sep = "\t",row.names = F,col.names = F)
