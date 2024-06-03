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
library(patchwork)
peak_preprocess_bin_level <- function(tissue,antibody,bin_size){
  tab = read.delim(paste0("data/samples/ATAC/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),skip=1)
  counts = tab[,c(7:10)]
  rownames(counts)= tab$Geneid
  colnames(counts) = c("young_1","old_1","young_2","old_2")
  group =c("young","old","young","old")
  y= DGEList(counts=counts,group=group)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$batch <- rep(c(rep("batch1", 2), rep("batch2", 2)), 1)
  y$samples$year <- group
  y$samples$year <- factor(y$samples$year,c("young","old"))
  y <- calcNormFactors(y)
  logCPMs <- cpm(y, log = TRUE)
  rv <- apply(logCPMs, 1, var)
  o <- order(rv, decreasing=TRUE)
  # top1000 <- head(o, 1000)
  # logCPM_top1000 <- logCPMs[top1000,]
  pca <- prcomp(t(logCPMs))
  to_plot <- data.frame(pca$x, y$samples)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  ggplot(to_plot, aes(x=PC1, y=PC2, color=batch, shape=year)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))
  # ggsave(paste0("result/",tissue,"/pca/",antibody,"_",bin_size,"_bins.png"),width = 5,height =5)
  
  batch <- factor(y$samples$batch)
  logCPMs_corrected <- limma::removeBatchEffect(logCPMs, batch = batch)
  # logCPM_corrected_top1000 <- logCPMs_corrected[top1000,]
  pca <- prcomp(t(logCPMs_corrected))
  to_plot <- data.frame(pca$x, y$samples)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  
  ggplot(to_plot, aes(x=PC1, y=PC2, color=batch, shape=year)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))
  # ggsave(paste0("result/",tissue,"/pca/",antibody,"_",bin_size,"_bins_pca_plot_after_remove_batch_effect.png"),width = 5,height =5)
  
  design <- model.matrix(~batch+year, y$samples)
  # y <- estimateDisp(y, design)
  
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = 3)
  tab<-tab[keep,]
  # fit <- glmQLFit(y, design)
  # fit  <- glmQLFTest(fit, coef = 3)
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
              "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC.old-young"=lrt$table$logFC)
  
  # txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  # # peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_macs2_mergecounts.bed"))
  # peak <- readPeakFile(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_merge-W",window_size,"-G",gap_size,"-E100.bed"))
  # peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
  #                          TxDb=txdb, annoDb="org.Mm.eg.db")
  # peakAnno<-unique(as.data.frame(peakAnno))
  # 
  # out$annotation <- peakAnno$annotation[which(peakAnno$V4 %in% out$Geneid)]
  # out$ens_id <- peakAnno$geneId[which(peakAnno$V4 %in% out$Geneid)]
  # out$symbol <- peakAnno$SYMBOL[which(peakAnno$V4 %in% out$Geneid)]
  out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                            ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
  out$Significant_bar <- "Stable"
  out$Significant_bar[which(out$`FDR.old-young` < 0.05 & (out$old_1/out$young_1 > 1.2) & (out$old_2/out$young_2 > 1.2))] <- "Up"
  out$Significant_bar[which(out$`FDR.old-young` < 0.05 & (out$old_1/out$young_1 < 0.8) & (out$old_2/out$young_2 < 0.8))] <- "Down"
  # colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
  # ggplot(
  #   # 数据、映射、颜色
  #   out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
  #   geom_point(aes(color = Significant), size=2) +
  #   scale_color_manual(values = colour[[nrow(as.data.frame(table(out$Significant)))]]) +
  #   # scale_color_manual(values = c("blue","grey")) +
  #   # 注释
  #   # geom_text_repel(
  #   #   data = subset(out,`FDR.FA-noFA` < 0.05 & abs(out$logFC) >= 1),
  #   #   aes(label = Geneid),
  #   #   size = 5,max.overlaps = 100,
  #   #   box.padding = unit(0.35, "lines"),
  #   #   point.padding = unit(0.3, "lines")) +
  #   # 辅助线
  #   geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  #   # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
  #   geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  #   # 坐标轴
  #   labs(x="log2(fold change)",
  #        y="-log10 (p-value)") +
  #   # xlim(-3,3)+
  #   # 图例
  #   theme_bw()+
  #   theme(text = element_text(size = 20))
  # ggsave(paste0("result/",tissue,"/diffpeaks/",antibody,"_merge-W",window_size,"-G",gap_size,"-E",e_value,"_volcano_plot_after_remove_batch_effect.png"),width = 10,height = 10)
  write.csv(out,paste0("data/samples/ATAC/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"),row.names = F)
  outup <- out[which(out$Significant=="Up"),]
  outdown <- out[which(out$Significant=="Down"),]
  write.table(outdown[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/ATAC/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect_down.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(outup[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/ATAC/",tissue,"/",antibody,"/bed/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect_up.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
antibodys <- c("ATAC")
# tissues <- c("stomach","skin")
# tissues <- c("aorta","tongue")
tissues <- c("bladder")
bin_size <-"10kb"
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody <- antibodys[j]
    peak_preprocess_bin_level(tissue,antibody,bin_size)
  }
}
bin_size <-"1kb"
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  for(j in c(1:length(antibodys))){
    antibody <- antibodys[j]
    peak_preprocess_bin_level(tissue,antibody,bin_size)
  }
}


bin_size <- "1kb"
tissues <- c("Hip","testis", "colon", "kidney", "lung", "spleen", "muscle", "cecum","bonemarrow","heart","thymus","aorta","tongue","skin","stomach","bladder")
p_list <- list()
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
# tissues <- c("heart","thymus")
for(k in c(1:length(tissues))){
  tissue <- tissues[k]
  out <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_1kb_bins_diff_after_remove_batch_effect.csv"))
  colour<- list(c("grey"),c("grey","red"),c("blue","grey","red"))
  quadrants<-c("first","second","third","fourth")
  quadrants_peaks <-list()
  for (i in c(1:length(quadrants))){
    quadrant <- quadrants[i]
    if(file.info(paste0("data/samples/ATAC/",tissue,"/ATAC/bed/H3K27me3_H3K9me3_all_significant_",bin_size,"_",quadrant,"_quadrant_unique_sort.bed"))$size >0){
      peaks<- read.delim(paste0("data/samples/ATAC/",tissue,"/ATAC/bed/H3K27me3_H3K9me3_all_significant_",bin_size,"_",quadrant,"_quadrant_unique_sort.bed"),header = F)
      quadrants_peaks[[i]] <- peaks$V4
      names(quadrants_peaks)[i] <- quadrant
    }
  }
  if(!is.null(quadrants_peaks[[1]])){
    p1<-ggplot() +
      geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",size = 0.5) +
      geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["first"]]),], mapping=aes(x = LogFC.old.young,y = -log10(FDR.old.young)),color = "#00b8a9",alpha=0.7,size = 0.5)+
      geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
      geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
      # 坐标轴
      labs(x="log2(fold change)",
           y="-log10 (FDR)") +
      # xlim(-3,3)+
      # 图例
      theme_bw()+
      theme(text = element_text(size = 10))
  }else{
    p1<-ggplot() +
      geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",size = 0.5) +
      geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
      geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
      # 坐标轴
      labs(x="log2(fold change)",
           y="-log10 (FDR)") +
      # xlim(-3,3)+
      # 图例
      theme_bw()+
      theme(text = element_text(size = 10))
  }
  if(!is.null(quadrants_peaks[[2]])){
    p2<-ggplot() +
      geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",size = 0.5) +
      geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["second"]]),], mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)),color ="#f6416c",alpha=0.7,size = 0.5)+
      geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
      geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
      # 坐标轴
      labs(x="log2(fold change)",
           y="-log10 (FDR)") +
      # xlim(-3,3)+
      # 图例
      theme_bw()+
      theme(text = element_text(size = 10))+
      annotate("text", x = min(out$LogFC.old.young), y = max(-log10(out$FDR.old.young)), label = nrow(out[which(out$Geneid %in%quadrants_peaks[["second"]] & out$LogFC.old.young<0),]), vjust = 5, hjust = 0,colour="blue",size=8)+
      annotate("text", x = max(out$LogFC.old.young), y = max(-log10(out$FDR.old.young)), label = nrow(out[which(out$Geneid %in%quadrants_peaks[["second"]] & out$LogFC.old.young>0),]), vjust = 5, hjust = 1.5,colour="red",size=8)
    
  }else{
    p2<-ggplot() +
      geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",size = 0.5) +
      geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
      geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
      # 坐标轴
      labs(x="log2(fold change)",
           y="-log10 (FDR)") +
      # xlim(-3,3)+
      # 图例
      theme_bw()+
      theme(text = element_text(size = 10))
  }
  if(!is.null(quadrants_peaks[[3]])){
  p3<-ggplot() +
    geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",size = 0.5) +
    geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["third"]]),], mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)),color ="#ffde7d",alpha=0.7,size = 0.5)+
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
    geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
    # 坐标轴
    labs(x="log2(fold change)",
         y="-log10 (FDR)") +
    # xlim(-3,3)+
    # 图例
    theme_bw()+
    theme(text = element_text(size = 10))
  }else{
    p3<-ggplot() +
      geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",size = 0.5) +
      geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
      geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
      # 坐标轴
      labs(x="log2(fold change)",
           y="-log10 (FDR)") +
      # xlim(-3,3)+
      # 图例
      theme_bw()+
      theme(text = element_text(size = 10))
  }
  if(!is.null(quadrants_peaks[[3]])){
  p4<-ggplot() +
    geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",size = 0.5) +
    geom_point(data=out[which(out$Geneid %in%quadrants_peaks[["fourth"]]),], mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)),color ="#48466d",alpha=0.7,size = 0.5)+
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
    geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
    # 坐标轴
    labs(x="log2(fold change)",
         y="-log10 (FDR)") +
    # xlim(-3,3)+
    # 图例
    theme_bw()+
    theme(text = element_text(size = 10))
  }else{
    p4<-ggplot() +
      geom_point(out, mapping=aes(x = LogFC.old.young, y = -log10(FDR.old.young)), color = "grey",size = 0.5) +
      geom_vline(xintercept=c(- 1,1),lty=4,col="black",lwd=0.8) +
      # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
      geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
      # 坐标轴
      labs(x="log2(fold change)",
           y="-log10 (FDR)") +
      # xlim(-3,3)+
      # 图例
      theme_bw()+
      theme(text = element_text(size = 10))
  }
  dir.create(paste0("result/",tissue,"/ATAC/"))
  if (tissue == "brain") tissue <- "FC"
  p2 <- p2+ ggtitle(tissue)
  p_list[[k]] <- p2
  # ggsave(paste0("result/",tissue,"/ATAC/diff_ATAC_with_H3K27me3_and_H3K9me3_relationship.png"),p,width = 4,height=4,type="cairo")
}
combined_plot <- plot_a_list(p_list, 4, 4)
ggsave(paste0("result/all/H3K27me3_H3K9me3/all_tissue_diff_ATAC_with_H3K27me3_and_H3K9me3_relationship.png"),combined_plot,width = 12,height=12,type="cairo")

bin_size <- "1kb"
tissues <- c("Hip","testis", "colon", "kidney", "lung", "spleen", "muscle", "cecum","bonemarrow","heart","thymus","aorta","tongue","skin","stomach","bladder")
peaks_anno_list_increase <- list()
peaks_anno_list_decrease <- list()
peaks_list_increase <- list()
peaks_list_decrease <- list()
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  out <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_1kb_bins_diff_after_remove_batch_effect.csv"))
  peaks<- read.delim(paste0("data/samples/ATAC/",tissue,"/ATAC/bed/H3K27me3_H3K9me3_all_significant_",bin_size,"_second_quadrant_unique_sort.bed"),header = F)
  out_increase <- out[which(out$Geneid %in% peaks$V4 & out$LogFC.old.young > 0),]
  peaks_list_increase[[i]] <- out_increase[,c("Chr","Start","End")]
  peaks_list_increase[[i]]$label <- paste0(peaks_list_increase[[i]]$Chr,"-",peaks_list_increase[[i]]$Start,"-",peaks_list_increase[[i]]$End)
  peaks_list_increase[[i]]$tissue <- tissue
  names(peaks_list_increase)[[i]] <- tissue
  peak_obj <- GRanges(seqnames = out_increase$Chr,   
                      ranges = IRanges(start = out_increase$Start, end = out_increase$End))
  peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  peaks_anno_list_increase[[i]] <- peak_anno
  names(peaks_anno_list_increase)[[i]] <- tissue
  
  out_decrease <- out[which(out$Geneid %in% peaks$V4 & out$LogFC.old.young < 0),]
  peaks_list_decrease[[i]] <- out_decrease[,c("Chr","Start","End")]
  peaks_list_decrease[[i]]$label <- paste0(peaks_list_decrease[[i]]$Chr,"-",peaks_list_decrease[[i]]$Start,"-",peaks_list_decrease[[i]]$End)
  peaks_list_decrease[[i]]$tissue <- tissue
  names(peaks_list_decrease)[[i]] <- tissue
  peak_obj <- GRanges(seqnames = out_decrease$Chr,   
                      ranges = IRanges(start = out_decrease$Start, end = out_decrease$End))
  peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  peaks_anno_list_decrease[[i]] <- peak_anno
  names(peaks_anno_list_decrease)[[i]] <- tissue
}

combined_increase_df <- do.call(rbind, peaks_list_increase) 
increase_merged_tissue <- combined_increase_df %>%   
  group_by(label) %>%   
  summarise(tissue = paste(tissue, collapse = "/")) %>%  
  ungroup() 
increase_label_counts <- combined_increase_df %>% count(label) 
table_combined_increase_df <- merge(increase_merged_tissue,increase_label_counts,by="label")
table_combined_increase_df <- table_combined_increase_df %>%  
  separate(label, into = c("Chr", "Start", "End"), sep = "-") 
table_combined_increase_df_select <- table_combined_increase_df[which(table_combined_increase_df$n>=8),]
table_combined_increase_df_select$Start <- as.numeric(table_combined_increase_df_select$Start)
table_combined_increase_df_select$End <- as.numeric(table_combined_increase_df_select$End)
peak_obj <- GRanges(seqnames = table_combined_increase_df_select$Chr,   
                    ranges = IRanges(start = table_combined_increase_df_select$Start, end = table_combined_increase_df_select$End))
peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Mm.eg.db")
peak_anno <- as.data.frame(peak_anno)
peak_anno <- peak_anno[,c(1:3,6,14,16,17)]

combined_decrease_df <- do.call(rbind, peaks_list_decrease) 
decrease_merged_tissue <- combined_decrease_df %>%   
  group_by(label) %>%   
  summarise(tissue = paste(tissue, collapse = "/")) %>%  
  ungroup() 
decrease_label_counts <- combined_decrease_df %>% count(label) 
table_combined_decrease_df <- merge(decrease_merged_tissue,decrease_label_counts,by="label")
table_combined_decrease_df <- table_combined_decrease_df %>%  
  separate(label, into = c("Chr", "Start", "End"), sep = "-") 
table_combined_decrease_df_select <- table_combined_decrease_df[which(table_combined_decrease_df$n>=3),]
table_combined_decrease_df_select$Start <- as.numeric(table_combined_decrease_df_select$Start)
table_combined_decrease_df_select$End <- as.numeric(table_combined_decrease_df_select$End)
peak_obj <- GRanges(seqnames = table_combined_decrease_df_select$Chr,   
                    ranges = IRanges(start = table_combined_decrease_df_select$Start, end = table_combined_decrease_df_select$End))
peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Mm.eg.db")
peak_anno <- as.data.frame(peak_anno)
peak_anno <- peak_anno[,c(1:3,6,14,16,17)]








peaks_anno_list_increase[[1]]@annoStat$Feature <- as.vector(peaks_anno_list_increase[[1]]@annoStat$Feature)
peaks_anno_list_increase[[1]]@annoStat$Feature <- factor(peaks_anno_list_increase[[1]]@annoStat$Feature,levels = c("Promoter (<=1kb)","Promoter (1-2kb)","Promoter (2-3kb)","5' UTR","3' UTR"  ,"1st Exon" ,"Other Exon","1st Intron"  , "Other Intron" ,"Distal Intergenic","Downstream (<=300)" ))
plotAnnoBar(peaks_anno_list_increase)

peaks_annolist_decrease[[1]]@annoStat$Feature <- as.vector(peaks_anno_list_decrease[[1]]@annoStat$Feature)
peaks_anno_list_decrease[[1]]@annoStat$Feature <- factor(peaks_anno_list_decrease[[1]]@annoStat$Feature,levels = c("Promoter (<=1kb)","Promoter (1-2kb)","Promoter (2-3kb)","5' UTR","3' UTR"  ,"1st Exon" ,"Other Exon","1st Intron"  , "Other Intron" ,"Distal Intergenic","Downstream (<=300)" ))
plotAnnoBar(peaks_anno_list_decrease)

for (i in c(16:length(tissues))){
  tissue <- tissues[i]
  up_peak_anno <- unique(as.data.frame(peaks_anno_list_increase[[i]]))
  genelist_up <- bitr(up_peak_anno$SYMBOL[which(str_detect(up_peak_anno$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  # genelist_up <- bitr(up_peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
  if(nrow(genelist_up_GO)>0){
    p <- barplot(genelist_up_GO)
    ggsave(paste0("result/all/H3K27me3_H3K9me3/ATAC_in_H3K27me3_up_H3K9me3_down/",tissue,"_increase_GO_promoter.png"),p,width=8,height = 10,type="cairo")
  }
  
  
  down_peak_anno <- unique(as.data.frame(peaks_anno_list_decrease[[i]]))
  # genelist_down <- bitr(down_peak_anno$SYMBOL[which(str_detect(down_peak_anno$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_down <- bitr(down_peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_down_GO <- enrichGO(genelist_down$ENTREZID,#GO富集分析
                               OrgDb = GO_database,
                               keyType = "ENTREZID",#设定读取的gene ID类型
                               ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                               pvalueCutoff = 0.05,#设定p值阈值
                               qvalueCutoff = 0.05,#设定q值阈值
                               readable = T)
  if(nrow(genelist_up_GO)>0){
    p <- barplot(genelist_down_GO)
    ggsave(paste0("result/all/H3K27me3_H3K9me3/ATAC_in_H3K27me3_up_H3K9me3_down/",tissue,"_decrease_GO.png"),p,width=8,height = 10,type="cairo")
  }
}
