rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
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
antibody = "H3K36me3"
tissue = "ileum"
# window_size = "1000"
# gap_size = "3000"
# e_value = "100"
# diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_merge-W1000-G3000-E100_diff_after_remove_batch_effect.csv"))
# diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_macs_narrowpeak_diff_after_remove_batch_effect.csv"))

bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    return ("10kb")
  }else{
    return("1kb")
  }  
}
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
antibodys <- c("H3K4me3","H3K4me1","H3K27ac")
tissues = c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus")
for(i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  for(j in c(1:length(tissues))){
    tissue <- tissues[j]
    out <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff_after_remove_batch_effect.csv"))
    out_up <- out[which(out$Significant_bar=="Up"),]
    out_down <- out[which(out$Significant_bar=="Down"),]
    if(nrow(out_up) > 0){
      up_peak <- GRanges(seqnames = out_up$Chr,   
                         ranges = IRanges(start = out_up$Start, end = out_up$End))
      up_peak_anno <- annotatePeak(up_peak, tssRegion=c(-3000, 3000),
                                   TxDb=txdb, annoDb="org.Mm.eg.db")
      up_peak_anno <- unique(as.data.frame(up_peak_anno))
      genelist_up <- bitr(up_peak_anno$SYMBOL[which(str_detect(up_peak_anno$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
      genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                                 OrgDb = GO_database,
                                 keyType = "ENTREZID",#设定读取的gene ID类型
                                 ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                                 pvalueCutoff = 0.05,#设定p值阈值
                                 qvalueCutoff = 0.05,#设定q值阈值
                                 readable = T)
      genelist_up_GO_result <- genelist_up_GO@result
      genelist_up_GO_result$tissue <- tissue
      if(j==1){
        up_anno_result <- genelist_up_GO_result
      }else{
        up_anno_result <- rbind(up_anno_result,genelist_up_GO_result)
      }
    }
    if(nrow(out_down) > 0){
      down_peak <- GRanges(seqnames = out_down$Chr,   
                           ranges = IRanges(start = out_down$Start, end = out_down$End))
      down_peak_anno <- annotatePeak(down_peak, tssRegion=c(-3000, 3000),
                                     TxDb=txdb, annoDb="org.Mm.eg.db")
      down_peak_anno <- unique(as.data.frame(down_peak_anno))
      genelist_down <- bitr(down_peak_anno$SYMBOL[which(str_detect(down_peak_anno$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
      genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                                  OrgDb = GO_database,
                                  keyType = "ENTREZID",#设定读取的gene ID类型
                                  ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                                  pvalueCutoff = 0.05,#设定p值阈值
                                  qvalueCutoff = 0.05,#设定q值阈值
                                  readable = T)
      genelist_down_GO_result <- genelist_down_GO@result
      genelist_down_GO_result$tissue <- tissue
      if(j==1){
        down_anno_result <- genelist_down_GO_result
      }else{
        down_anno_result <- rbind(down_anno_result,genelist_down_GO_result)
      }
    }
  }
  up_anno_result_table<-as.data.frame(table(up_anno_result$ID[which(up_anno_result$p.adjust<0.05)]))
  down_anno_result_table<-as.data.frame(table(down_anno_result$ID[which(down_anno_result$p.adjust<0.05)]))

  up_anno_result_filter <- up_anno_result[which(up_anno_result$ID %in% up_anno_result_table$Var1[which(up_anno_result_table$Freq>=4)]),c("Description","p.adjust","tissue")]
  up_anno_result_plot <- reshape2::melt(up_anno_result_filter)
  up_anno_result_plot$value <- -log10(up_anno_result_plot$value)
  ggplot(up_anno_result_plot, aes(x=tissue, y=Description, fill=value)) +
    geom_tile()+
    scale_fill_gradient2(low="#a8d8ea", high="#990033",mid = "#ffffd2")+
    theme(axis.text.x = element_text(angle=90, vjust=0.5),axis.text = element_text(size = 12)) +
    labs(title = paste(antibody,"up"),x=NULL, y=NULL, fill="-log10(p.adjust)")+theme_bw()
  ggsave(paste0("result/all/diff/",antibody,"/bins_promoter_increase_GO_heatmap.png"),width = 10,height = 15,type="cairo")
  down_anno_result_filter <- down_anno_result[which(down_anno_result$ID %in% down_anno_result_table$Var1[which(down_anno_result_table$Freq>=4)]),c("Description","p.adjust","tissue")]
  down_anno_result_plot <- reshape2::melt(down_anno_result_filter)
  down_anno_result_plot$value <- -log10(down_anno_result_plot$value)
  ggplot(down_anno_result_plot, aes(x=tissue, y=Description, fill=value)) +
    geom_tile()+
    scale_fill_gradient2(low="#a8d8ea", high="#990033",mid = "#ffffd2")+
    theme(axis.text.x = element_text(angle=90, vjust=0.5),axis.text = element_text(size = 12)) +
    labs(title = paste(antibody,"down"),x=NULL, y=NULL, fill="-log10(p.adjust)")+theme_bw()
  ggsave(paste0("result/all/diff/",antibody,"/bins_promoter_decrease_GO_heatmap.png"),width = 10,height = 15,type="cairo")
  
}


















tissues <- c("spleen","testis","colon","kidney","lung","liver","muscle","Hip","cecum","bonemarrow","brain","ileum","heart","thymus")
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  diff <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff_after_remove_batch_effect.csv"))
  diff_up <- diff[which(diff$Significant_bar=="Up"),]
  peak_obj <- GRanges(seqnames = diff_up$Chr,   
                      ranges = IRanges(start = diff_up$Start, end = diff_up$End))
  peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  peak_anno <- unique(as.data.frame(peak_anno))
  # genelist_up <-bitr(diff$symbol[which(diff$Significant_bar=="Up" & str_detect(diff$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  # genelist_up <-bitr(diff$symbol[which(diff$Significant_bar=="Up")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  # genelist_up <-bitr(diff$symbol[which(diff$FDR.old.young<0.05 & diff$LogFC.old.young > 0.5)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_up <- bitr(peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_up_GO <-enrichGO( genelist_up$ENTREZID,#GO富集分析
                             OrgDb = GO_database,
                             keyType = "ENTREZID",#设定读取的gene ID类型
                             ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                             pvalueCutoff = 0.05,#设定p值阈值
                             qvalueCutoff = 0.05,#设定q值阈值
                             readable = T)
  if(nrow(genelist_up_GO) > 0){
    barplot(genelist_up_GO,font.size=15,label_format = 100)
    dir.create(paste0("result/",tissue,"/",antibody))
    ggsave(paste0("result/",tissue,"/",antibody,"/increase_Go.png"),width = 8,height = 7,type="cairo")
    GO_up <- genelist_up_GO@result 
    write.csv(GO_up,paste0("result/",tissue,"/",antibody,"/GO_up.csv"))
  }
  # genelist_down <-bitr(diff$symbol[which(diff$Significant_bar=="Down" & str_detect(diff$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  # genelist_down <-bitr(diff$symbol[which(diff$FDR.old.young<0.05 & diff$LogFC.old.young < -0.5)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  diff_down <- diff[which(diff$Significant_bar=="Down"),]
  peak_obj <- GRanges(seqnames = diff_down$Chr,   
                      ranges = IRanges(start = diff_down$Start, end = diff_down$End))
  peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  peak_anno <- unique(as.data.frame(peak_anno))
  genelist_down <- bitr(peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_down_GO <-enrichGO( genelist_down$ENTREZID,#GO富集分析
                               OrgDb = GO_database,
                               keyType = "ENTREZID",#设定读取的gene ID类型
                               ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                               pvalueCutoff = 0.05,#设定p值阈值
                               qvalueCutoff = 0.05,#设定q值阈值
                               readable = T)
  if(nrow(genelist_down_GO) > 0){
    barplot(genelist_down_GO,font.size=15,label_format = 100)
    dir.create(paste0("result/",tissue,"/",antibody))
    ggsave(paste0("result/",tissue,"/",antibody,"/decrease_Go.png"),width = 8,height = 7,type="cairo")
    GO_down <- genelist_down_GO@result 
    write.csv(GO_up,paste0("result/",tissue,"/",antibody,"/GO_down.csv"))
  }
}

antibody <- "H3K4me3"
GO_up <- list()
GO_down <- list()
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  if (file.exists(paste0("result/",tissue,"/",antibody,"/GO_up.csv"))) {  
    diff <- read.csv(paste0("result/",tissue,"/",antibody,"/GO_up.csv"))
    GO_up[[i]]<-diff
    names(GO_up)[[i]] <- tissue
  }
  if (file.exists(paste0("result/",tissue,"/",antibody,"/GO_down.csv"))) {  
    diff <- read.csv(paste0("result/",tissue,"/",antibody,"/GO_down.csv"))
    GO_down[[i]]<-diff
    names(GO_down)[[i]] <- tissue
  }
}
GO_up <- Filter(Negate(is.null), GO_up)
GO_down <- Filter(Negate(is.null), GO_down)

GO_up_data <- Reduce(function(x, y) merge(x, y, by = "ID", all = F), lapply(GO_up, function(x) x[, c("ID", "Description","p.adjust")]))
GO_up_description <- GO_up_data[,seq(2, ncol(GO_up_data), by = 2)]
description <- GO_up_description[,1]
for (i in c(2:ncol(GO_up_description))){
  description <- coalesce(description,GO_up_description[,i])
}
GO_up_data <- GO_up_data[,-seq(2, ncol(GO_up_data), by = 2)]
GO_up_data$description <- description
colnames(GO_up_data)[2:(ncol(GO_up_data)-1)] <- names(GO_up)
GO_up_data[,2:(ncol(GO_up_data)-1)] <- -log10(GO_up_data[,2:(ncol(GO_up_data)-1)])
# GO_up_data$count <- rowSums(GO_up_data[,2:(ncol(GO_up_data)-1)] > -log10(0.05), na.rm = TRUE)  
# GO_up_data <- GO_up_data[GO_up_data$count >= 4, ] 
GO_up_data_plot <- reshape2::melt(GO_up_data)
rownames(GO_up_data) <- GO_up_data$ID
GO_up_data <- GO_up_data[,-1]
set.seed(1)
kmeans_centers <- 4
kmeans <- kmeans(GO_up_data[,-c(ncol(GO_up_data))],centers=kmeans_centers)
cluster<-as.data.frame(kmeans$cluster)
cluster$rownames<-rownames(cluster)
cluster<-cluster[order(cluster$`kmeans$cluster`),]
GO_up_data_plot$ID <- factor(GO_up_data_plot$ID,levels=cluster$rownames)
colnames(cluster)[2]<-"ID"
GO_up_data_plot <- merge(GO_up_data_plot,cluster,by="ID")
colnames(GO_up_data_plot)[5] <- "kmeans"
GO_up_data_plot$kmeans <- factor(GO_up_data_plot$kmeans,levels=seq(kmeans_centers, 1, -1)  )
ggplot(GO_up_data_plot, aes(x=variable, y=ID, fill=value)) +
  geom_tile(color = "black") +
  facet_grid(rows = vars(kmeans),space="free_y",scales="free_y")+
  scale_fill_gradient2(limits = c(0, -log(0.001)),low="#3490de", 
                       mid="#f6f6f6", high="#f73859",  midpoint = -log10(0.05),na.value="black",oob = scales::squish)  +
  theme(axis.text.x = element_text(angle=90, vjust=0.5),axis.text = element_text(size = 12)) +
  scale_y_discrete(breaks=NULL)+
  labs(title = paste(antibody,"up"),x=NULL, y=NULL, fill=NULL)
kmeans4 <- GO_up_data$description[which(rownames(GO_up_data) %in% cluster$ID[which(cluster$`kmeans$cluster` == 4)])]
kmeans3 <- GO_up_data$description[which(rownames(GO_up_data) %in% cluster$ID[which(cluster$`kmeans$cluster` == 3)])]
kmeans2 <- GO_up_data$description[which(rownames(GO_up_data) %in% cluster$ID[which(cluster$`kmeans$cluster` == 2)])]
kmeans1 <- GO_up_data$description[which(rownames(GO_up_data) %in% cluster$ID[which(cluster$`kmeans$cluster` == 1)])]

GO_down_data <- Reduce(function(x, y) merge(x, y, by = "ID", all = F), lapply(GO_down, function(x) x[, c("ID", "Description","p.adjust")]))
GO_down_description <- GO_down_data[,seq(2, ncol(GO_down_data), by = 2)]
description <- GO_down_description[,1]
for (i in c(2:ncol(GO_down_description))){
  description <- coalesce(description,GO_down_description[,i])
}
GO_down_data <- GO_down_data[,-seq(2, ncol(GO_down_data), by = 2)]
GO_down_data$description <- description
colnames(GO_down_data)[2:(ncol(GO_down_data)-1)] <- names(GO_down)
GO_down_data[,2:(ncol(GO_down_data)-1)] <- -log10(GO_down_data[,2:(ncol(GO_down_data)-1)])
# GO_down_data$count <- rowSums(GO_down_data[,2:(ncol(GO_down_data)-1)] > -log10(0.05), na.rm = TRUE)  
# GO_down_data <- GO_down_data[GO_down_data$count >= 4, ] 
GO_down_data_plot <- reshape2::melt(GO_down_data)
rownames(GO_down_data) <- GO_down_data$ID
GO_down_data <- GO_down_data[,-1]
set.seed(1)
kmeans_centers <- 4
kmeans <- kmeans(GO_down_data[,-c(ncol(GO_down_data))],centers=kmeans_centers)
cluster<-as.data.frame(kmeans$cluster)
cluster$rownames<-rownames(cluster)
cluster<-cluster[order(cluster$`kmeans$cluster`),]
GO_down_data_plot$ID <- factor(GO_down_data_plot$ID,levels=cluster$rownames)
colnames(cluster)[2]<-"ID"
GO_down_data_plot <- merge(GO_down_data_plot,cluster,by="ID")
colnames(GO_down_data_plot)[5] <- "kmeans"
GO_down_data_plot$kmeans <- factor(GO_down_data_plot$kmeans,levels=seq(kmeans_centers, 1, -1)  )
ggplot(GO_down_data_plot, aes(x=variable, y=ID, fill=value)) +
  geom_tile(color = "black") +
  facet_grid(rows = vars(kmeans),space="free_y",scales="free_y")+
  scale_fill_gradient2(limits = c(0, -log(0.001)),low="#3490de", 
                       mid="#f6f6f6", high="#f73859",  midpoint = -log10(0.05),na.value="black",oob = scales::squish)  +
  theme(axis.text.x = element_text(angle=90, vjust=0.5),axis.text = element_text(size = 12)) +
  scale_y_discrete(breaks=NULL)+
  labs(title = paste(antibody,"down"),x=NULL, y=NULL, fill=NULL)
kmeans4 <- GO_down_data$description[which(rownames(GO_down_data) %in% cluster$ID[which(cluster$`kmeans$cluster` == 4)])]
kmeans3 <- GO_down_data$description[which(rownames(GO_down_data) %in% cluster$ID[which(cluster$`kmeans$cluster` == 3)])]
kmeans2 <- GO_down_data$description[which(rownames(GO_down_data) %in% cluster$ID[which(cluster$`kmeans$cluster` == 2)])]
kmeans1 <- GO_down_data$description[which(rownames(GO_down_data) %in% cluster$ID[which(cluster$`kmeans$cluster` == 1)])]



