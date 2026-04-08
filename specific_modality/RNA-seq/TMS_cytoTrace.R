rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(CytoTRACE2)
tab = read.csv("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/GSE132040/GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv")
tab <- tab[-c(54353:54357),]
rownames(tab) <- tab[,1]
tab <- tab[,-1]
new_col_name <- sub("\\..*", "", colnames(tab)) 
colnames(tab) <- new_col_name
cytotrace2_result <- cytotrace2(tab,species = "mouse",slot_type="counts",seed=1,ncores=10)



table = read.delim("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/GSE132040/GSE132040_MACA_Bulk_metadata.csv",sep = ',')
table$tissue <- sub("_\\d.*$", "", table$source.name)
cytotrace2_result <- merge(cytotrace2_result,table[,c(1,3,5,7,14)],by.x="row.names",by.y="Sample.name")
cytotrace2_result <- cytotrace2_result[which(cytotrace2_result$characteristics..sex=="m"),]
to_plot <- cytotrace2_result[,c("Row.names","CytoTRACE2_Score","tissue","characteristics..age")]
average_scores <-aggregate(CytoTRACE2_Score ~ tissue, data = to_plot, FUN = mean)
average_scores <- average_scores[order(average_scores$CytoTRACE2_Score),]
to_plot$tissue <- factor(to_plot$tissue,levels = average_scores$tissue)
color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sort(unique(as.character(to_plot$tissue))))
ggplot(to_plot,aes(x=tissue,y=CytoTRACE2_Score,color = tissue))+    
  geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7)+
  scale_color_manual(values = color)+
  ggtitle(paste0("Stemness Score"))+ylim(0,1)+
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("")+labs(fill = "", color = "") +ylab("CytoTrace2 stemness score")
