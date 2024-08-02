rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)

# 设置tissue变量
if (length(args) > 0) {
  tissue <- args[1]
} else {
  tissue <- "default value"  # 如果没有提供参数，使用默认值
}
# tissue <- "brain"
H3K27me3 <- read.delim(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/H3K27me3_percentage_intersect_H3K27me3_08_H3K9me3_08_merge.bed"),
                       header = F)
H3K9me3 <- read.delim(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/H3K9me3_percentage_intersect_H3K27me3_08_H3K9me3_08_merge.bed"),
                      header = F)

H3K27me3$label <- paste(H3K27me3$V1,H3K27me3$V2,H3K27me3$V3,sep="-")
H3K9me3$label <- paste(H3K9me3$V1,H3K9me3$V2,H3K9me3$V3,sep="-")
data <- merge (H3K27me3[,7:8],H3K9me3[,7:8],by="label")
colnames(data)[2:3]<-c("H3K27me3","H3K9me3")

p<- ggplot(data, aes(H3K9me3, H3K27me3))+ 
  stat_density2d(aes(fill=..density..), geom="raster", contour=FALSE)+
  scale_fill_gradient(low = "#ffffd2", high = "#f38181", limits = c(0.0001, 0.0008), oob = scales::squish)+
  theme_minimal()+
  theme(text = element_text(size = 15))+
  labs(x = "% H3K9me3 overlap", y = "% H3K27me3 overlap")
ggsave(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/peak_level/distribution.png"),p,width = 7,height = 6,type="cairo")
colocalized <- data[which(data$H3K27me3>75),]
embedded <- data[which(data$H3K27me3<=75 & data$H3K9me3 >=75),]
overlap <- data[which(data$H3K27me3>=75 & data$H3K9me3 >=75),]

write.table(H3K27me3[which(H3K27me3$label %in% colocalized$label),c(1:3)],
            paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/peak_level/colocalized_tmp.bed"),sep = "\t",quote = F,row.names = F,col.names = F)
write.table(H3K9me3[which(H3K9me3$label %in% embedded$label),c(1:3)],
            paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/peak_level/embedded.bed"),sep="\t",quote = F,row.names = F,col.names = F)

write.table(H3K9me3[which(H3K9me3$label %in% overlap$label),c(1:3)],
            paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/peak_level/overlap.bed"),sep="\t",quote = F,row.names = F,col.names = F)
