rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
WGBS_R1 <- read.delim("data/raw_data/20240804_WGBS/XX315/XX315_S1_L002_R1/bed/XX315_S1_L002_R1_CpG.bdg",skip = 1,header = F)
WGBS_R2 <- read.delim("data/raw_data/20240804_WGBS/XX315/XX315_S1_L002_R2/bed/XX315_S1_L002_R2_CTCpG.bdg",skip = 1,header = F)
head(WGBS_R1)
data_list <- list(WGBS_R1=WGBS_R1,
             WGBS_R2=WGBS_R2)
for(i in c(1:2)){
  data_list[[i]]$depth <- data_list[[i]]$V5
  data_list[[i]]$percent <- round(data_list[[i]]$V4/data_list[[i]]$V5,2)*100
}
for(i in c(1:length(data_list))){
  data_list[[i]] <- data_list[[i]][which(data_list[[i]]$depth > 15),]
  data_list[[i]]$label <- paste0(data_list[[i]]$V1,"-",data_list[[i]]$V2,"-",data_list[[i]]$V3)
}
i=1
j=2
df <- merge(data_list[[i]][,c("label","percent")],data_list[[j]][,c("label","percent")],by="label")
colnames(df) <- c("label",names(data_list)[[i]],names(data_list)[[j]])
par(cex.lab = 1.5, cex.axis = 1.2)  
smoothScatter(df[,2] ~ df[,3],xlab = names(data_list)[[j]],ylab = names(data_list)[[i]])
abline(a = 0, b = 1, col = "red", lty = 2)  
mtext(paste0(names(data_list)[[i]]," = ", round(mean(df[,2]),2)), side = 3, line = -2, adj = 0.05,  cex = 1.2)  
mtext(paste0(names(data_list)[[j]]," = ", round(mean(df[,3]),2)), side = 3, line = -3.5, adj = 0.05,  cex = 1.2)  
mtext(paste0("r = ", round(cor(df[,2],df[,3]),2)), side = 1, line = -1.5, adj = 0.9,  cex = 1.2)  
