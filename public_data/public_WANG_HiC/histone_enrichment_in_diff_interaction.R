rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(reshape2)
library(data.table)
bed <- read.table("data/public_data/WANG_cellular_aging_HiC/raw_matrix/DS1_200000_abs.bed")
bed <- bed[which(bed$V1 %in% paste0("chr",c(1:22,"X"))),]
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
# SRR560 <- read.table("data/public_data/cellular_aging_GSE133292/Chip_seq/peaks/macs_narrowpeak/SRR13274560_10kb.txt")
# SRR561 <- read.table("data/public_data/cellular_aging_GSE133292/Chip_seq/peaks/macs_narrowpeak/SRR13274561_10kb.txt")
ENCFF265JIG <-  read.table("data/public_data/NIH_Roadmap_BJ/peaks/macs_narrowpeak/ENCFF265JIG_10kb.txt")
# counts <- cbind(SRR560,SRR561[,"V4"])
# counts$average <- (counts[,4]+counts[,5])/2
counts <- ENCFF265JIG
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


re <-  read.table(paste0("data/public_data/WANG_cellular_aging_HiC/differential_analysis/G_DS_200000.FDR"))
re$Significant <- "Stable"
re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 < 0)] <- "Down"
re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 > 0)] <- "Up"
# re_sig <- re[which(re$Significant!="Stable"),]
re_sig <- re
re_sig <- re_sig[which(abs(re_sig$V2 - re_sig$V3)>4),]
re_sig$V1 <- paste0("chr",re_sig$V1)
re_sig$V1[which(re_sig$V1=="chr23")] <- "chrX"
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
# tab <- read.table("data/public_data/cellular_aging_GSE133292/Chip_seq/H3K9me3_G_20kb_bins.counts",header = T)
# counts = tab[,c(7:ncol(tab))]
# rownames(counts)= tab$Geneid
# counts <- cpm(counts)
# counts_average <- rowSums(counts)/2
# counts <- cbind(tab[,c(1:6)],counts_average)
# counts$label <- paste0(counts$Chr,"-",counts$Start,"-",counts$End)
# 
# re_sig <- merge(re_sig,counts[,c(7:8)],by.x="region1",by.y="label")
# colnames(re_sig)[which(colnames(re_sig)=="counts_average")] <- "region1_average"
# re_sig <- merge(re_sig,counts[,c(7:8)],by.x="region2",by.y="label")
# colnames(re_sig)[which(colnames(re_sig)=="counts_average")] <- "region2_average"
# SRR560 <- read.table("data/public_data/cellular_aging_GSE133292/Chip_seq/peaks/macs_narrowpeak/SRR13274560_20kb_board.txt")
# SRR561 <- read.table("data/public_data/cellular_aging_GSE133292/Chip_seq/peaks/macs_narrowpeak/SRR13274561_20kb_board.txt")
# SRR569 <- read.table("data/public_data/cellular_aging_GSE133292/Chip_seq/peaks/macs_narrowpeak/SRR13274569_20kb_board.txt")
# counts <- cbind(SRR560,SRR561[,"V4"])
# counts$average <- (counts[,4]+counts[,5])
# counts$V1 <- paste0("chr",counts$V1)
# counts$V1[which(counts$V1=="chr23")] <- "chrX"
# counts$V2 <- counts$V2-1
# counts <- data.frame(label=paste(counts$V1,counts$V2,counts$V3,sep = "-"), average = counts$average)
# re_sig <- merge(re_sig,counts,by.x="region1",by.y="label")
# colnames(re_sig)[which(colnames(re_sig)=="average")] <- "region1_average"
# re_sig <- merge(re_sig,counts,by.x="region2",by.y="label")
# colnames(re_sig)[which(colnames(re_sig)=="average")] <- "region2_average"
# re_sig$Significant <- factor(re_sig$Significant,levels=c("Up","Down","Stable"))
ggplot()+
  geom_point(data=re_sig[which(re_sig$Significant == "Down"),], mapping=aes(region1_average, region2_average),color="blue",size=0.1) +
  geom_point(data=re_sig[which(re_sig$Significant == "Up"),], mapping=aes(region1_average, region2_average),color="red",size=0.1) +
  theme_bw()+xlab("Enrichment in anchor1")+ylab("Enrichment in anchor2")+
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black") +  
  geom_hline(yintercept = 2.5, linetype = "dashed", color = "black") +  
  ggtitle("H3K9me3 in G")+ ylim(0,5) + xlim(0,5)+
  theme(text = element_text(size = 14))  

type1 <- re_sig[which(re_sig$Significant != "Stable" & re_sig$region1_average>1.5 & re_sig$region2_average>1.5),]
type2 <- re_sig[which(re_sig$Significant != "Stable" & re_sig$region1_average<1.5 & re_sig$region2_average>1.5),]
type3 <- re_sig[which(re_sig$Significant != "Stable" & re_sig$region1_average<1.5 & re_sig$region2_average<1.5),]
type1 <- as.data.frame(table(type1$Significant))
type1$condition <- "type1"
type1$percent <- type1$Freq/sum(type1$Freq)*100
type2 <- as.data.frame(table(type2$Significant))
type2$condition <- "type2"
type2$percent <- type2$Freq/sum(type2$Freq)*100
type3 <- as.data.frame(table(type3$Significant))
type3$condition <- "type3"
type3$percent <- type3$Freq/sum(type3$Freq)*100
to_plot <- rbind(type1,type2,type3)
color <- setNames(c("red","blue"),c("Up","Down"))
ggplot(to_plot, aes(x = condition, y = percent, fill = Var1)) +  
  geom_bar(stat = 'identity',color="white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  ylab("Proportion")

color <- setNames(c("red","blue"),c("Up","Down"))

re_sig_plot <- re_sig[which(re_sig$Significant != "Stable"),]
re_sig_plot <- sample_n(re_sig_plot,50000)
ggplot(data=re_sig_plot, aes(region1_average, region2_average,color=Significant))+
  geom_point(size=0.05) +
  scale_color_manual(values = color)+
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "black") +  
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "black") + 
  theme_bw()+xlab("Enrichment in anchor1")+ylab("Enrichment in anchor2")+
  ggtitle("H3K9me3 in G") + ylim(0,5) + xlim(0,5)

# ggplot(re_sig[which(re_sig$Significant != "Stable"),], aes(x = region1_average)) +   
#   geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +  
#   labs(title = "Sum 10kb in 200kb", x = "Value", y = "Frequency") +  
#   theme_minimal()+ xlim(0,10)  
# 
# ggplot(SRR560, aes(x = V4)) +   
#   geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +  
#   labs(title = "G rep1 10kb enrichment", x = "Value", y = "Frequency") +  
#   theme_minimal() + xlim(0,1)
# 
# ggplot(SRR561, aes(x = V4)) +   
#   geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +  
#   labs(title = "G rep2 10kb enrichment", x = "Value", y = "Frequency") +  
#   theme_minimal() + xlim(0,1)

