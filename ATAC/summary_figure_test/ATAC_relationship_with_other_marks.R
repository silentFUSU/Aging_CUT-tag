rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
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
#####WGBS
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
summary <- list()
for(tissue in tissues){
  diff_peak <- list()
  df <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set_DAR/",tissue,"_DARs.txt"))
  peaks <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set/",tissue,".bed"))
  rownames(peaks) <- paste0(peaks$V1,":",peaks$V2,"-",peaks$V3)
  df <- merge(df,peaks,by="row.names")
  colnames(df)[1] <- "label"
  diff_peak[["increase"]] <- df[which(df$logFC > 0 & df$FDR < 0.05),]
  diff_peak[["decrease"]] <- df[which(df$logFC < 0 & df$FDR < 0.05),]
  search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
  t_search_table <- search_table[which(search_table$tissue==tissue),]
  tissue_summary <- list(increase=data.frame(),decrease=data.frame())
  for(condition in c("increase","decrease")){
    t_diff_peak <- diff_peak[[condition]]
    if(nrow(t_diff_peak) >= 100){
      t_diff_peak <- as.data.table(t_diff_peak)
      setDT(t_diff_peak)
      setkey(t_diff_peak,V1,V2,V3)
      for(sample in t_search_table$sample_name){
        df_WGBS <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
        setDT(df_WGBS)
        setkey(df_WGBS,V1,V2,V3)  
        overlaps <- foverlaps(df_WGBS, t_diff_peak, type = "any", nomatch = 0L)  
        
        result <- overlaps[, .(V4_sum = sum(V4), V5_sum = sum(V5)), by = label]
        result <- as.data.frame(result)
        result$methylation <- result$V4_sum/result$V5_sum
        result <- result[,c("label","methylation")]
        colnames(result)[2] <- sample
        if(nrow(tissue_summary[[condition]])==0){
          tissue_summary[[condition]] <- result
        }else{
          tissue_summary[[condition]] <- merge(tissue_summary[[condition]],result,by="label")
        }
      }
      t_tissue_summary <- tissue_summary[[condition]]
      young_summary <- t_tissue_summary[,c("label",t_search_table$sample_name[which(t_search_table$age=="3M")])]
      old_summary <- t_tissue_summary[,c("label",t_search_table$sample_name[which(t_search_table$age=="24M")])]
      young_summary$young_methylation <- rowMeans(young_summary[,-1])
      old_summary$old_methylation <- rowMeans(old_summary[,-1])
      t_tissue_summary <- merge(young_summary,old_summary,by="label")
      t_tissue_summary$delta <- t_tissue_summary$old_methylation - t_tissue_summary$young_methylation
      t_tissue_summary <- t_tissue_summary[,c("label","delta")]
      tissue_summary[[condition]] <- merge(tissue_summary[[condition]],t_tissue_summary,by="label")
    }
  }
  summary[[tissue]] <- tissue_summary
}
summary <- readRDS("data/samples/ATAC/all/ATAC/WGBS_change_in_ATAC_diff_peak_from_LMJ.rds")
conditions <- c("increase","decrease")
to_plot <- data.frame()
for(tissue in tissues){
  for(condition in conditions){
    df <- summary[[tissue]][[condition]]
    if(nrow(df) > 0){
      t_summary <- data.frame(tissue=tissue_label_change(tissue),delta=median(df$delta),condition=condition)
    }else{
      t_summary <- data.frame(tissue=tissue_label_change(tissue),delta=NA,condition=condition)
    }
    to_plot <- rbind(to_plot,t_summary)
  }  
}
tissues_order <- c("Ovary","iWAT","Mammary Gland","BAT","Thymus","Uterus","Bone Marrow","Spleen","Liver","Skin","Lung","Aorta","Bladder","Cerebellum","Muscle","Hippocampus","Pancreas","Jejunum","Kidney","Cortex","Colon","Tongue","Heart","Cecum","Stomach","Ileum","Testis")
to_plot$tissue <- factor(to_plot$tissue,levels = rev(tissues_order))
to_plot$delta[which(to_plot$delta > 0.2)] <- 0.2
to_plot$delta[which(to_plot$delta < -0.2)] <- -0.2

ggplot(to_plot, aes(x = condition, y = tissue, fill = delta)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",limits = c(-0.2, 0.2), midpoint = 0) +
  theme_minimal() +
  ggtitle(paste0("DNA methylation delta in DARs"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### histone modification
## correlation
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")

for(tissue in tissues){
  ATAC_search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
  ATAC_df <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_1kb_bins.counts"),header = T)
  rownames(ATAC_df) <- ATAC_df$Geneid
  
  ATAC_peaks <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set/",tissue,".bed"))
  setDT(ATAC_peaks)
  setkey(ATAC_peaks,V1,V2,V3)
  
  ATAC_bin <- ATAC_df[,c(1:4)]
  ATAC_bin$Start <- ATAC_bin$Start+1
  ATAC_bin <- as.data.table(ATAC_bin)
  setDT(ATAC_bin)
  setkey(ATAC_bin,Chr,Start,End)
  
  overlaps <- as.data.frame(foverlaps(ATAC_bin,ATAC_peaks, type = "any", nomatch = 0L))
  
  ATAC_df <- ATAC_df[which(rownames(ATAC_df) %in% overlaps$Geneid),]
  
  ATAC_summary <- read.table(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_1kb_bins.counts.summary"),header = T)
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+|HM[0-9]+).*"
  colnames(ATAC_df)[c(7:ncol(ATAC_df))] <- gsub(pattern, "\\1",colnames(ATAC_df)[c(7:ncol(ATAC_df))])
  colnames(ATAC_summary) <- gsub(pattern,"\\1",colnames(ATAC_summary))
  ATAC_search_table <- ATAC_search_table[which(ATAC_search_table$sample_name %in% colnames(ATAC_df)),]
  
  rownames(ATAC_df) <- ATAC_df$Geneid
  ATAC_df <- ATAC_df[,c("Chr","Start","End","Length",ATAC_search_table$sample_name)]
  ATAC_counts <- ATAC_df[,ATAC_search_table$sample_name]
  ATAC_summary <- ATAC_summary[-2,ATAC_search_table$sample_name]
  total_reads <- colSums(ATAC_summary)
  length <- as.numeric(ATAC_df$Length)
  rpkm <- sweep(ATAC_counts,2,total_reads,"/")
  rpkm <- sweep(rpkm,1,length,"/") * 1000000000
  
  rpkm_young <- rpkm[,ATAC_search_table$sample_name[which(ATAC_search_table$age=="3m")]]
  rpkm_old <- rpkm[,ATAC_search_table$sample_name[which(ATAC_search_table$age=="24m")]]
  rpkm_young$mean_young <- rowMeans(rpkm_young)
  rpkm_old$mean_old <- rowMeans(rpkm_old)
  rpkm_mean_summary <- merge(rpkm_young[,"mean_young",drop=F],rpkm_old[,"mean_old",drop=F],by="row.names")
  
  rpkm_mean_summary$log2FC <- log2(rpkm_mean_summary$mean_old/rpkm_mean_summary$mean_young)
  colnames(rpkm_mean_summary)[which(colnames(rpkm_mean_summary)=="log2FC")] <-"ATAC"
  colnames(rpkm_mean_summary)[1] <- "Geneid"
  
  ATAC_rpkm_summary <- rpkm_mean_summary[,c("Geneid","ATAC")]
  
  for(antibody in antibodys){
    CUTTAG_search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
    df <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins.counts"),header = T)
    rownames(df) <- df$Geneid
    summary <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins.counts.summary"),header = T)
    pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+|HM[0-9]+).*"
    colnames(df)[c(7:ncol(df))] <- gsub(pattern, "\\1",colnames(df)[c(7:ncol(df))])
    colnames(summary) <- gsub(pattern,"\\1",colnames(summary))
    CUTTAG_search_table <- CUTTAG_search_table[which(CUTTAG_search_table$sample_name %in% colnames(df)),]
    counts <- df[,CUTTAG_search_table$sample_name]
    summary <- summary[-2,CUTTAG_search_table$sample_name]
    
    total_reads <- colSums(summary)
    length <- as.numeric(df$Length)
    rpkm <- sweep(counts,2,total_reads,"/")
    rpkm <- sweep(rpkm,1,length,"/") * 1000000000
    rpkm_young <- rpkm[,CUTTAG_search_table$sample_name[which(CUTTAG_search_table$age=="3m")]]
    rpkm_old <- rpkm[,CUTTAG_search_table$sample_name[which(CUTTAG_search_table$age=="24m")]]
    rpkm_young$mean_young <- rowMeans(rpkm_young)
    rpkm_old$mean_old <- rowMeans(rpkm_old)
    rpkm_mean_summary <- merge(rpkm_young[,"mean_young",drop=F],rpkm_old[,"mean_old",drop=F],by="row.names")
    rpkm_mean_summary <- rpkm_mean_summary[which(rpkm_mean_summary$Row.names %in% ATAC_rpkm_summary$Geneid),]

    rpkm_mean_summary$log2FC <- log2(rpkm_mean_summary$mean_old/rpkm_mean_summary$mean_young)
    rpkm_mean_summary <- rpkm_mean_summary[!is.na(rpkm_mean_summary$log2FC) & !is.infinite(rpkm_mean_summary$log2FC), ]
    colnames(rpkm_mean_summary)[which(colnames(rpkm_mean_summary)=="log2FC")] <- antibody
    colnames(rpkm_mean_summary)[1] <- "Geneid"
    cor_df <- merge(ATAC_rpkm_summary,rpkm_mean_summary[,c("Geneid",antibody)],by="Geneid")
    cor_test <- cor.test(cor_df[,2],cor_df[,3])
    t_cor_summary <- data.frame(tissue=tissue_label_change(tissue),antibody=antibody,cor=as.numeric(cor_test$estimate["cor"]))
    cor_summary <- rbind(cor_summary,t_cor_summary)
  }
}

## median logFC
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
tissues <-  c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
              "thymus","skin","bladder","bonemarrow","Hip","heart",
              "muscle","jejunum","uterus","ovary","liver","tongue",
              "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
conditions <- c("increase","decrease")
antibody_logFC_summary <- list()
for(antibody in antibodys){
  logFC_summary <- list()
  for(tissue in tissues){
    diff_peak <- list()
    df <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set_DAR/",tissue,"_DARs.txt"))
    peaks <- read.table(paste0("data/samples/ATAC/ATAC_peak_from_LMJ/tissue_peak_set/",tissue,".bed"))
    rownames(peaks) <- paste0(peaks$V1,":",peaks$V2,"-",peaks$V3)
    df <- merge(df,peaks,by="row.names")
    colnames(df)[1] <- "label"
    diff_peak[["increase"]] <- df[which(df$logFC > 0 & df$FDR < 0.05),]
    diff_peak[["decrease"]] <- df[which(df$logFC < 0 & df$FDR < 0.05),]
    search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
    t_search_table <- search_table[which(search_table$tissue==tissue),]
    tissue_summary <- list(increase=data.frame(),decrease=data.frame())
    for(condition in c("increase","decrease")){
      t_diff_peak <- diff_peak[[condition]]
      if(nrow(t_diff_peak) >= 100){
        t_diff_peak <- as.data.table(t_diff_peak)
        setDT(t_diff_peak)
        setkey(t_diff_peak,V1,V2,V3)
        
        CUTTAG_search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
        df <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins.counts"),header = T)
        rownames(df) <- df$Geneid
        df_bin <- df[,c(1:4)]
        df_bin <- as.data.table(df_bin)
        setDT(df_bin)
        setkey(df_bin,Chr,Start,End)
        overlaps <- as.data.frame(foverlaps(df_bin,t_diff_peak, type = "any", nomatch = 0L))
        
        summary <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins.counts.summary"),header = T)
        pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+|HM[0-9]+).*"
        colnames(df)[c(7:ncol(df))] <- gsub(pattern, "\\1",colnames(df)[c(7:ncol(df))])
        colnames(summary) <- gsub(pattern,"\\1",colnames(summary))
        CUTTAG_search_table <- CUTTAG_search_table[which(CUTTAG_search_table$sample_name %in% colnames(df)),]
        counts <- df[,CUTTAG_search_table$sample_name]
        summary <- summary[-2,CUTTAG_search_table$sample_name]
        
        total_reads <- colSums(summary)
        length <- as.numeric(df$Length)
        rpkm <- sweep(counts,2,total_reads,"/")
        rpkm <- sweep(rpkm,1,length,"/") * 1000000000
        rpkm_young <- rpkm[,CUTTAG_search_table$sample_name[which(CUTTAG_search_table$age=="3m")]]
        rpkm_old <- rpkm[,CUTTAG_search_table$sample_name[which(CUTTAG_search_table$age=="24m")]]
        rpkm_young$mean_young <- rowMeans(rpkm_young)
        rpkm_old$mean_old <- rowMeans(rpkm_old)
        rpkm_mean_summary <- merge(rpkm_young[,"mean_young",drop=F],rpkm_old[,"mean_old",drop=F],by="row.names")
        
        rpkm_mean_summary <- rpkm_mean_summary[which(rpkm_mean_summary$Row.names %in% overlaps$Geneid),]
        
        rpkm_mean_summary$log2FC <- log2(rpkm_mean_summary$mean_old/rpkm_mean_summary$mean_young)
        tissue_summary[[condition]] <- rpkm_mean_summary
      }
    }
    logFC_summary[[tissue]] <- tissue_summary
  }
  antibody_logFC_summary[[antibody]] <- logFC_summary
}

conditions <- c("increase","decrease")
for(antibody in antibodys){
  logFC_summary <- antibody_logFC_summary[[antibody]]
  to_plot <- data.frame()
  for(tissue in tissues){
    for(condition in conditions){
      df <- logFC_summary[[tissue]][[condition]]
      if(nrow(df) > 0){
        t_logFC_summary <- data.frame(tissue=tissue_label_change(tissue),log2FC=median(df$log2FC,na.rm = T),condition=condition)
      }else{
        t_logFC_summary <- data.frame(tissue=tissue_label_change(tissue),log2FC=NA,condition=condition)
      }
      to_plot <- rbind(to_plot,t_logFC_summary)
    }  
  }
  tissues_order <- c("Ovary","iWAT","Mammary Gland","BAT","Thymus","Uterus","Bone Marrow","Spleen","Liver","Skin","Lung","Aorta","Bladder","Cerebellum","Muscle","Hippocampus","Pancreas","Jejunum","Kidney","Cortex","Colon","Tongue","Heart","Cecum","Stomach","Ileum","Testis")
  to_plot$tissue <- factor(to_plot$tissue,levels = rev(tissues_order))
  to_plot$log2FC[which(to_plot$log2FC > 1)] <- 1
  to_plot$log2FC[which(to_plot$log2FC < -1)] <- -1
  ggplot(to_plot, aes(x = condition, y = tissue, fill = log2FC)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",limits = c(-1, 1), midpoint = 0) +
    theme_minimal() +
    ggtitle(paste0(antibody," log2FC in DARs"))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
