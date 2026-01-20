rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
library(ggsignif)
library(data.table)
library(dplyr)
library(GenomeInfoDb)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library("GenomicRanges")
library(genomation)
library(ChIPseeker)
library(clusterProfiler)
library(enrichplot)
tissue_label_change <- function(tissue){
  if(tissue=="brain"){
    tissue_label <- "Cortex"
  }else if(tissue == "Hip"){
    tissue_label <- "Hippocampus"
  }else if(tissue == "CB"){
    tissue_label <- "Cerebellum"
  }else{
    tissue_label <- str_to_title(tissue)
    if(tissue_label == "Bonemarrow"){
      tissue_label <- "Bone Marrow"
    }else if(tissue_label == "Bat"){
      tissue_label <- "BAT"
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
tissues <- c("ileum","jejunum","colon","cecum") 
PRC2_LMR_list <- list()
LMR_condition_summary <- data.frame()
for(tissue in tissues){
  hmr <- read.table(paste0("data/samples/WGBS/",tissue,"/hmr/all_samples_hmr.bed"))
  if(tissue %in% c("mammarygland","ovary","uterus")){
    hmr <- hmr[which(hmr$V1 %in% paste0("chr",c(1:19,"X"))),]
  }else{
    hmr <- hmr[which(hmr$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  }
  
  blacklist <- read.table("~/ref_data/mm10-blacklist.v2.bed",sep = "\t")
  blacklist <- as.data.table(blacklist)
  setDT(blacklist)
  setkey(blacklist,V1,V2,V3)
  
  hmr_regions <- as.data.table(hmr)
  setDT(hmr_regions)
  setkey(hmr_regions,V1,V2,V3)
  overlaps <- foverlaps(hmr_regions, blacklist, type = "any", nomatch = 0L)  
  
  hmr <- hmr[which(!hmr$V4 %in% overlaps$i.V4),]
  
  EZH2 <- read.table("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_input/EZH2_peaks.narrowPeak")
  EZH2 <- EZH2[,c(1:3,9)]
  colnames(EZH2)[4] <- "EZH2"
  
  SUZ12 <- read.table("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_input/SUZ12_peaks.narrowPeak")
  SUZ12 <- SUZ12[,c(1:3,9)]
  colnames(SUZ12)[4] <- "SUZ12"
  
  hmr_regions <- as.data.table(hmr)
  setDT(hmr_regions)
  setkey(hmr_regions,V1,V2,V3)
  
  EZH2 <- as.data.table(EZH2)
  setDT(EZH2)
  setkey(EZH2,V1,V2,V3)
  
  SUZ12 <- as.data.table(SUZ12)
  setDT(SUZ12)
  setkey(SUZ12,V1,V2,V3)
  
  overlaps <- foverlaps(EZH2,hmr_regions, type = "any", nomatch = 0L)  
  overlaps <- as.data.table(overlaps[,c("V1","V2","V3","V4","EZH2")])
  setDT(overlaps)
  setkey(overlaps,V1,V2,V3)
  
  overlaps_SUZ12 <- foverlaps(overlaps, SUZ12, type = "any", nomatch = 0L)
  overlaps_SUZ12$SUZ12 <- as.numeric(overlaps_SUZ12$SUZ12)
  overlaps_SUZ12$EZH2 <- as.numeric(overlaps_SUZ12$EZH2)
  overlaps_SUZ12$score <- (overlaps_SUZ12$SUZ12 + overlaps_SUZ12$EZH2)/2
  score <- overlaps_SUZ12 %>%
    group_by(V4) %>%
    summarize(
      score = max(score, na.rm = TRUE)
    )
  score <- as.data.frame(score)
  score <- score[order(score$score,decreasing = T),]
  score_top <- score[1:min(1000,nrow(score)),]
  
  regions <- as.data.table(hmr[which(hmr$V4 %in%  score_top$V4),c(1:4)])
  
  setDT(regions)
  setkey(regions,V1,V2,V3)
  
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  summary <- data.frame()
  for(sample in search_table$sample_name){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    setDT(df)
    setkey(df,V1,V2,V3)  
    overlaps <- foverlaps(df, regions, type = "any", nomatch = 0L)  
    result <- overlaps[, .(V4_sum = sum(i.V4), V5_sum = sum(V5)), by = V4]
    result <- as.data.frame(result)
    result$methylation <- result$V4_sum/result$V5_sum*100
    result <- result[which(result$V5_sum >=20),]
    result <- result[,c("V4","methylation")]
    
    colnames(result)[2] <- sample 
    if(nrow(summary)==0){
      summary <- result
    }else{
      summary <- merge(summary,result,by="V4")
    }
  }
  young_summary <- summary[,c("V4",search_table$sample_name[which(search_table$age=="3M")])]
  young_summary$young <- rowMeans(young_summary[,-1])
  
  old_summary <- summary[,c("V4",search_table$sample_name[which(search_table$age=="24M")])]
  old_summary$old <- rowMeans(old_summary[,-1])
  
  summary <- merge(young_summary[,c("V4","young")],old_summary[,c("V4","old")],by="V4")
  summary$delta <- summary$old - summary$young
  summary$condition <- "Up"
  summary$condition[which(summary$delta < 0 )] <- "Down"
  
  PRC2_LMR <- hmr[which(hmr$V4 %in% summary$V4[which(summary$condition=="Up")]),c(1:3)]
  PRC2_LMR_list[[tissue]] <- PRC2_LMR
  
  table <- as.data.frame(table(summary$condition))
  t_LMR_condition_summary <- data.frame(tissue=tissue_label_change(tissue),Up=table$Freq[which(table$Var1=="Up")],Down=table$Freq[which(table$Var1=="Down")])
  LMR_condition_summary <- rbind(LMR_condition_summary,t_LMR_condition_summary)
  }

gr_ileum <- GRanges(seqnames = PRC2_LMR_list$ileum$V1,
                    ranges = IRanges(start = PRC2_LMR_list$ileum$V2,
                                     end = PRC2_LMR_list$ileum$V3))

gr_jejunum <- GRanges(seqnames = PRC2_LMR_list$jejunum$V1,
                      ranges = IRanges(start = PRC2_LMR_list$jejunum$V2,
                                       end = PRC2_LMR_list$jejunum$V3))

gr_colon <- GRanges(seqnames = PRC2_LMR_list$colon$V1,
                    ranges = IRanges(start = PRC2_LMR_list$colon$V2,
                                     end = PRC2_LMR_list$colon$V3))

gr_cecum <- GRanges(seqnames = PRC2_LMR_list$cecum$V1,
                    ranges = IRanges(start = PRC2_LMR_list$cecum$V2,
                                     end = PRC2_LMR_list$cecum$V3))

overlaps_ij <- findOverlaps(gr_ileum, gr_jejunum)
overlaps_ijc <- findOverlaps(gr_ileum[queryHits(overlaps_ij)], gr_colon)
overlaps_all <- findOverlaps(gr_ileum[queryHits(overlaps_ijc)], gr_cecum)

gr_list <- list(ileum=gr_ileum,jejunum=gr_jejunum,colon=gr_colon,cecum=gr_cecum)

overlap_summary <- data.frame()
for(tissue in tissues){
  overlapping_regions <- as.data.frame(gr_list[[tissue]][queryHits(overlaps_all)])
  t_overlap_summary <- data.frame(tissue=tissue_label_change(tissue),overlap=nrow(overlapping_regions),all=nrow(PRC2_LMR_list[[tissue]]))
  t_overlap_summary$percentage <- t_overlap_summary$overlap/t_overlap_summary$all *100
  overlap_summary <- rbind(overlap_summary,t_overlap_summary)
  }

overlapping_regions <- as.data.frame(gr_list[["ileum"]][queryHits(overlaps_all)])
write.table(overlapping_regions[,1:3],"data/samples/WGBS/all/LMR_intestine_PRC2_hypermethylation.bed",append = F,quote = F,sep = "\t",row.names = F,col.names = F)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
gr <- GRanges(seqnames = overlapping_regions$seqnames,
              ranges = IRanges(start = overlapping_regions$start, end = overlapping_regions$end))
peakAnno <- annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db")
peakAnno_df <- as.data.frame(peakAnno)
plotAnnoBar(peakAnno)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'

genelist <- bitr(unique(peakAnno_df$SYMBOL),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
genelist_GO <- enrichGO( genelist$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
barplot(genelist_GO,label_format = 50,showCategory = 20)

genelist_GO <- pairwise_termsim(genelist_GO)
emapplot(genelist_GO, showCategory = 50) 

### background
EZH2 <- read.table("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_input/EZH2_peaks.narrowPeak")
EZH2 <- EZH2[,c(1:4,9)]
colnames(EZH2)[5] <- "EZH2"

SUZ12 <- read.table("data/public_data/EZH2_SUZ12_E14/peaks/macs_narrowpeak_input/SUZ12_peaks.narrowPeak")
SUZ12 <- SUZ12[,c(1:4,9)]
colnames(SUZ12)[5] <- "SUZ12"

EZH2 <- as.data.table(EZH2)
setDT(EZH2)
setkey(EZH2,V1,V2,V3)

SUZ12 <- as.data.table(SUZ12)
setDT(SUZ12)
setkey(SUZ12,V1,V2,V3)
overlaps <- foverlaps(EZH2,SUZ12, type = "any", nomatch = 0L)  
score <- overlaps %>%
  group_by(V1,V2,V3,V4,SUZ12) %>%
  summarize(
    EZH2 = max(EZH2, na.rm = TRUE)
  )
score <- as.data.frame(score)
score$score <- (score$SUZ12+score$EZH2)/2
score <-score[order(score$score,decreasing = T),]
score <- score[c(1:1000),]
# score <- score
bg_gr <- GRanges(seqnames = score$V1,
              ranges = IRanges(start = score$V2, end = score$V3))
peakAnno <- annotatePeak(bg_gr, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db")
peakAnno_df <- as.data.frame(peakAnno)
plotAnnoBar(peakAnno)
genelist_bg <- bitr(unique(peakAnno_df$SYMBOL),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)

genelist_GO <- enrichGO( genelist$ENTREZID,
                         univers=genelist_bg$ENTREZID,
                         OrgDb = GO_database,
                         keyType = "ENTREZID",#设定读取的gene ID类型
                         ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                         pvalueCutoff = 0.05,#设定p值阈值
                         qvalueCutoff = 0.05,#设定q值阈值
                         readable = T)
barplot(genelist_GO,label_format = 50,showCategory = 20)
