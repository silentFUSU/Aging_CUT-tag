rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
convert2bedpe <- function(re_sig,resolution){
  df <- re_sig[,c("region1","region2","Significant")]
  df <- df %>%  
    separate(region1, into = c("chr1", "x1", "x2"), sep = "-", convert = TRUE)  
  df <- df %>%  
    separate(region2, into = c("chr2", "y1", "y2"), sep = "-", convert = TRUE)  
  write.table(df[which(df$Significant=="Up"),c(1:6)],paste0("data/samples/HiC/",tissue,"/differential_analysis/Wang_output_",resolution,"_increase.bedpe"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(df[which(df$Significant=="Down"),c(1:6)],paste0("data/samples/HiC/",tissue,"/differential_analysis/Wang_output_",resolution,"_decrease.bedpe"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  
}

tissue <- "brain"
antibody <- "H3K9me3"
search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff.csv")
search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
search_table <- search_table[which(search_table$age=="3m"),]
young1 <- read.table(paste0("data/samples/HiC/",tissue,"/with_histone/",antibody,"/",search_table$sample_name[1],"_10kb.txt"))
young2 <- read.table(paste0("data/samples/HiC/",tissue,"/with_histone/",antibody,"/",search_table$sample_name[2],"_10kb.txt"))

resolution <- "200000"
HiC_search_table <- read.csv("data/samples/all/HiC_search_table.csv")
HiC_search_table <- HiC_search_table[which(HiC_search_table$tissue==tissue),]
bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",HiC_search_table$sample_name[1],"_",resolution,"_abs.bed"))
bed <- bed[which(bed$V1 %in% paste0("chr",c(1:19,"X"))),]
bed$V5 <- "NA"
count <- 1
bed[1,"V5"] <- count
for(i in c(2:nrow(bed))){
  if(bed[i,"V1"] != bed[i-1,"V1"]){
    count <- 1
    bed[i,"V5"] <- count
  }else{
    count <- count+1
    bed[i,"V5"] <- count
  }
}

counts <- cbind(young1,young2[,"V4"])
counts$average <- (counts[,4]+counts[,5])/2
counts$average <- counts[,4]
counts$V1 <- paste0("chr",counts$V1)
counts$V1[which(counts$V1=="chr23")] <- "chrX"
counts <- counts[,c("V1","V2","V3","average")]
bed$V2 <- bed$V2+1
bed <- as.data.table(bed)
setDT(bed)
setkey(bed,V1,V2,V3)
counts <- as.data.table(counts)
setDT(counts)
setkey(counts,V1,V2,V3)
overlaps <- foverlaps(counts, bed, type = "any", nomatch = 0L)  
overlaps$V2 <- overlaps$V2-1
bed <- overlaps
bed$region <- paste(bed$V1,bed$V2,bed$V3,sep = "-")
bed$label <- paste(bed$V1,bed$V5,sep = "-")
bed <- as.data.frame(bed)
bed <- bed %>%  
  group_by(label) %>%  
  mutate(enrichment = sum(average)) %>%  
  ungroup() %>%  
  select(region, label, enrichment) %>%  
  distinct(label, .keep_all = TRUE)  

re <-  read.table(paste0("data/samples/HiC/",tissue,"/differential_analysis/",tissue,"_",resolution,".FDR"))
re$Significant <- "Stable"
re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 < 0)] <- "Down"
re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 > 0)] <- "Up"
re_sig <- re[which(re$Significant!="Stable"),]

# re_sig <- re
re_sig <- re_sig[which(abs(re_sig$V2 - re_sig$V3)>4),]
re_sig$V1 <- paste0("chr",re_sig$V1)
re_sig$V1[which(re_sig$V1=="chr20")] <- "chrX"
re_sig_reverse <- data.frame(V1=re_sig$V1,V2=re_sig$V3, V3=re_sig$V2, V4=re_sig$V4,V5=re_sig$V5,V6=re_sig$V6,Significant=re_sig$Significant)
re_sig <- rbind(re_sig,re_sig_reverse)

re_sig$label1 <- paste(re_sig$V1,re_sig$V2,sep = "-")
re_sig$label2 <- paste(re_sig$V1,re_sig$V3,sep = "-")
re_sig <- merge(re_sig,bed,by.x="label1",by.y="label")
colnames(re_sig)[which(colnames(re_sig)=="region")] <- "region1"
re_sig <- merge(re_sig,bed,by.x="label2",by.y="label")
colnames(re_sig)[which(colnames(re_sig)=="region")] <- "region2"
colnames(re_sig)[which(colnames(re_sig)=="enrichment.x")] <- "region1_average"
colnames(re_sig)[which(colnames(re_sig)=="enrichment.y")] <- "region2_average"

# ggplot()+
#   geom_point(data=re_sig[which(re_sig$Significant == "Down"),], mapping=aes(region1_average, region2_average),color="blue",size=0.1) +
#   geom_point(data=re_sig[which(re_sig$Significant == "Up"),], mapping=aes(region1_average, region2_average),color="red",size=0.1) +
#   theme_bw()+xlab("Enrichment in anchor1")+ylab("Enrichment in anchor2")+
#   geom_vline(xintercept = 2.5, linetype = "dashed", color = "black") +  
#   geom_hline(yintercept = 2.5, linetype = "dashed", color = "black") +  
#   ggtitle("H3K9me3 in G")+ ylim(0,10) + xlim(0,10)+
#   theme(text = element_text(size = 14))  
color <- setNames(c("red","blue"),c("Up","Down"))

re_sig_plot <- sample_n(re_sig,min(50000,nrow(re_sig)))
ggplot(data=re_sig_plot, aes(region1_average, region2_average,color=Significant))+
  geom_point(size=0.05) +
  scale_color_manual(values = color)+
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "black") +  
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "black") + 
  theme_bw()+xlab("Enrichment in anchor1")+ylab("Enrichment in anchor2")+
  ggtitle(paste0(antibody," in Young")) + ylim(0,10) + xlim(0,10)



