rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(ChIPseeker)
library(EnsDb.Mmusculus.v79)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(maditr)
state <- c("E13")
target_state <- "E6"
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin") 
dir.create(paste0("result/all/ChromHMM/15_until_skin/state_transfer/",paste0(state, collapse="-"),"_to_",target_state))
dir.create(paste0("result/all/ChromHMM/15_until_skin/state_transfer/",paste0(state, collapse="-"),"_to_",target_state,"/bed"))
peaks_anno_list <-list()
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
peaks_list <- list()
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  young1_bed <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_young1_15_segments_1k.bed"),header = F)
  young2_bed <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_young2_15_segments_1k.bed"),header = F)
  old1_bed <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_old1_15_segments_1k.bed"),header = F)
  old2_bed <- read.delim(paste0("result/all/ChromHMM/15_until_skin/split_1k/",tissue,"_old2_15_segments_1k.bed"),header = F)
  young1_bed <- young1_bed[which(young1_bed$V4%in%state),]
  young2_bed <- young2_bed[which(young2_bed$V4%in%state),]
  young_bed <- intersect(young1_bed,young2_bed)
  old1_bed <- old1_bed[which(old1_bed$V4==target_state),]
  old2_bed <- old2_bed[which(old2_bed$V4==target_state),]
  old_bed <- intersect(old1_bed,old2_bed)
  young_bed$label <- paste0(young_bed$V1,"-",young_bed$V2,"-",young_bed$V3)
  old_bed$label <- paste0(old_bed$V1,"-",old_bed$V2,"-",old_bed$V3)
  region_interest <- merge(young_bed,old_bed,by="label")
  region_interest <- young_bed[which(young_bed$label %in% old_bed$label),c(1:3)]
  write.table(region_interest,
              file=paste0(paste0("result/all/ChromHMM/15_until_skin/state_transfer/",paste0(state, collapse="-"),"_to_",target_state,"/bed/",tissue,".bed")),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  peaks_list[[i]] <- region_interest
  peaks_list[[i]]$label <- paste0(peaks_list[[i]]$V1,"-",peaks_list[[i]]$V2,"-",peaks_list[[i]]$V3)
  peaks_list[[i]]$tissue <- tissue
  names(peaks_list)[[i]] <- tissue
  peak_obj <- GRanges(seqnames = region_interest$V1,   
                      ranges = IRanges(start = region_interest$V2, end = region_interest$V3))
  peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  peaks_anno_list[[i]] <- peak_anno
  names(peaks_anno_list)[[i]] <- tissue
  }
names(peaks_anno_list)[[1]] <- "FC"
peaks_anno_list$FC@annoStat$Feature <- as.vector(peaks_anno_list$FC@annoStat$Feature)
peaks_anno_list$FC@annoStat$Feature <- factor(peaks_anno_list$FC@annoStat$Feature,levels = c("Promoter (<=1kb)","Promoter (1-2kb)","Promoter (2-3kb)","5' UTR","3' UTR"  ,"1st Exon" ,"Other Exon","1st Intron"  , "Other Intron" ,"Distal Intergenic","Downstream (<=300)" ))
plotAnnoBar(peaks_anno_list)

combined_peak <- do.call(rbind, peaks_list) 
merged_tissue <- combined_peak %>%   
  group_by(label) %>%   
  summarise(tissue = paste(tissue, collapse = "/")) %>%  
  ungroup() 
label_counts <- combined_peak %>% count(label) 
table_combined_peak <- merge(merged_tissue,label_counts,by="label")
gtf_data <- read.delim("/storage/zhangyanxiaoLab/share/gtf/mm10.gencode.vM25.annotation.gtf", header = FALSE, skip = 5)  
gtf_data <- gtf_data[which(gtf_data$V3=="gene"),]
gtf_data <- separate(gtf_data, V9, into = paste0("V9_", 1:6), sep = "; ") 
gtf_data <- gtf_data[,c(1:5,11)]
gtf_data$length <- gtf_data$V5-gtf_data$V4+1
gtf_data$V9_3 <- gsub("gene_name ", "", gtf_data$V9_3) 
colnames(gtf_data)[6] <- "Symbol"
test <- vector()
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  if (tissue == "brain") tissue <- "FC"
  peak_anno <- as.data.frame(peaks_anno_list[[tissue]])
  genelist <- bitr(peak_anno$SYMBO,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  gtf_data_tissue <- gtf_data[which(gtf_data$Symbol %in% genelist$SYMBOL),]
  gtf_data_other <- gtf_data[which(! gtf_data$Symbol %in% genelist$SYMBOL),]
  # gtf_data_other_random <- gtf_data_other[sample(nrow(gtf_data_other), nrow(gtf_data_tissue)),]
  gtf_data_tissue$condition <- "cryptic transcription"
  # gtf_data_other_random$condition <- "random"
  gtf_data_other$condition <- "other all"
  gtf_data_tissue$tissue <- tissue
  gtf_data_other$tissue <- tissue
  t_plot <- rbind(gtf_data_tissue,gtf_data_other)
  if (i == 1) {
    plot<-t_plot
  }else{
    plot<-rbind(plot,t_plot)
  }
  t_test <- t.test(gtf_data_tissue$length, gtf_data_other$length) 
  test <- c(test,t_test$p.value)
}

tissues[1] <-"FC"
plot$tissue <- factor(plot$tissue,levels = tissues)
ggplot(plot,aes(x=tissue,y=log10(length),fill =condition))+
  geom_boxplot()+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("log10(gene length)")+
  xlab("")+labs(fill = "", color = "")+theme(axis.text.x = element_text(angle = 45, hjust = 1))  


t.test(gtf_data_tissue$length, gtf_data_other$length) 
ggplot(plot,aes(x=condition,y=length,fill =condition))+
  geom_boxplot()+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("gene length")+
  xlab("")+labs(fill = "", color = "")

genelist_GO <-enrichGO( genelist$ENTREZID,#GO富集分析
                           OrgDb = GO_database,
                           keyType = "ENTREZID",#设定读取的gene ID类型
                           ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                           pvalueCutoff = 0.05,#设定p值阈值
                           qvalueCutoff = 0.05,#设定q值阈值
                           readable = T)
barplot(genelist_GO,showCategory = 10,label_format = 50)
