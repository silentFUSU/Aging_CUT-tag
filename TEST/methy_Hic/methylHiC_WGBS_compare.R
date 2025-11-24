rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)

juicer_101 <- read.delim("data/raw_data/20240725_methylHiC/lambda_test_all_data/WJH_Mousebrain_C1_BSseq/WJH_Mousebrain_C1_BSseq_101/aligned/merged_dedup_sort_CpG.bedGraph",skip = 1,header=F)
juicer_102 <- read.delim("data/raw_data/20240725_methylHiC/lambda_test_all_data/WJH_Mousebrain_C1_BSseq/WJH_Mousebrain_C1_BSseq_102/aligned/merged_dedup_sort_CpG.bedGraph",skip = 1,header=F)
bhmem_101 <- read.delim("data/raw_data/20240725_methylHiC/bhmem/WJH_Mousebrain_C1_BSseq_101/vcf/WJH_Mousebrain_C1_BSseq_101.calmd.nodup.cpg.raw.sort.CG.6plus2.bed",skip = 1,header=F)
bhmem_102 <- read.delim("data/raw_data/20240725_methylHiC/bhmem/WJH_Mousebrain_C1_BSseq_102/vcf/WJH_Mousebrain_C1_BSseq_102.calmd.nodup.cpg.raw.sort.CG.6plus2.bed",skip = 1,header=F)
WGBS_101 <- read.delim("data/raw_data/20240804_WGBS/XX315/XX315_S1_L002/bed/XX315_S1_L002_CpG.bdg",skip = 1,header=F)
WGBS_102 <- read.delim("data/raw_data/20240804_WGBS/XX316/XX316_S2_L002/bed/XX316_S2_L002_CpG.bdg",skip = 1,header=F)

data_list <- list(juicer_101=juicer_101,juicer_102=juicer_102,
                  bhmem_101=bhmem_101,bhmem_102=bhmem_102,
                  WGBS_101=WGBS_101,WGBS_102=WGBS_102)


for(i in c(1:2)){
  data_list[[i]]$depth <- data_list[[i]]$V5+data_list[[i]]$V6
  data_list[[i]]$percent <- round(data_list[[i]]$V5/data_list[[i]]$depth,2)*100
}

for(i in c(3:4)){
  data_list[[i]]$depth <- data_list[[i]]$V8
  data_list[[i]]$percent <- data_list[[i]]$V7
}

for(i in c(5:6)){
  data_list[[i]]$depth <- data_list[[i]]$V5
  data_list[[i]]$percent <- round(data_list[[i]]$V4/data_list[[i]]$V5,2)*100
}

for(i in c(1:length(data_list))){
  data_list[[i]] <- data_list[[i]][which(data_list[[i]]$depth > 15),]
  data_list[[i]]$label <- paste0(data_list[[i]]$V1,"-",data_list[[i]]$V2,"-",data_list[[i]]$V3)
}


i=4
j=6
df <- merge(data_list[[i]][,c("label","percent")],data_list[[j]][,c("label","percent")],by="label")
colnames(df) <- c("label",names(data_list)[[i]],names(data_list)[[j]])
par(cex.lab = 1.5, cex.axis = 1.2)  
smoothScatter(df[,2] ~ df[,3],xlab = names(data_list)[[j]],ylab = names(data_list)[[i]])
abline(a = 0, b = 1, col = "red", lty = 2)  
mtext(paste0(names(data_list)[[i]]," = ", round(mean(df[,2]),2)), side = 3, line = -2, adj = 0.05,  cex = 1.2)  
mtext(paste0(names(data_list)[[j]]," = ", round(mean(df[,3]),2)), side = 3, line = -3.5, adj = 0.05,  cex = 1.2)  
mtext(paste0("r = ", round(cor(df[,2],df[,3]),2)), side = 1, line = -1.5, adj = 0.9,  cex = 1.2)  

i=1
j=3
k=5
heatmap_df <- merge(data_list[[i]][,c("label","percent")],data_list[[j]][,c("label","percent")],by="label")
heatmap_df <- merge(heatmap_df,data_list[[k]][,c("label","percent")],by="label")
colnames(heatmap_df) <- c("label",names(data_list)[[i]],names(data_list)[[j]],names(data_list)[[k]])
pheatmap::pheatmap(heatmap_df[2:4],cluster_cols = F,show_rownames = F)
save.image("data/raw_data/20240725_methylHiC/WGBS_methylHIC-compare")
