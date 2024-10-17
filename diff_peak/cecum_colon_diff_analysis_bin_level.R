rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(ggrepel)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me3","H3K4me1")
PCA_analysis <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
  cecum <- read.delim(paste0("data/samples/cecum/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),skip = 1)
  colon <- read.delim(paste0("data/samples/colon/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),skip = 1)
  tab <- merge(cecum,colon[,c(1,c(7:ncol(colon)))],by="Geneid")
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  print(colnames(counts))
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
  search_table <- search_table[order(search_table$sample_name),]
  age <- search_table$age
  mouse_ID <- search_table$mouse_ID
  tissue <- search_table$tissue
  age[which(age=="3m")] <- "young"
  age[which(age=="24m")] <- "old"
  # colnames(counts) <- paste0(colnames(counts),"-",age,"-",mouse_ID)
  batch <- c("batch1","batch2","batch2","batch3","batch1","batch3","batch4","batch4")
  y <- DGEList(counts=counts)  
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  logCPM <- cpm(y,log = T)
  logCPM_corrected <- limma::removeBatchEffect(logCPM,batch = batch)
  pca <- prcomp(t(logCPM))
  to_plot <- data.frame(pca$x)
  to_plot$sample_name <- rownames(to_plot)
  to_plot <- merge(to_plot,search_table[,-2],by="sample_name")
  to_plot$sample_name <- paste0(to_plot$sample_name,"-",to_plot$age,"-",to_plot$mouse_ID)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  to_plot$age <- factor(to_plot$age,levels = c("3m","24m"))
  ggplot(to_plot, aes(x=PC1, y=PC2, color=tissue,shape=age)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
    geom_text_repel(  
      data = to_plot,  
      aes(x = PC1, y = PC2, label = sample_name, color = tissue),  
      size = 5,  
      box.padding = unit(0.35, "lines"),  
      point.padding = unit(0.3, "lines")  
    )+
    ggtitle(antibody)
}
p_list <- list()
for(i in c(1:length(antibodys))){
  p_list[[i]] <- PCA_analysis(antibodys[i]) 
}
combined_plot <- plot_a_list(p_list,2,3)
ggsave("result/intestine/colon_cecum_PCA.png",combined_plot,width = 18,height = 10,type="cairo")

diff_analysis <- function(antibody){
  p_list <- list(colon=list(),cecum=list())
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
  cecum <- read.delim(paste0("data/samples/cecum/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),skip = 1)
  colon <- read.delim(paste0("data/samples/colon/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),skip = 1)
  tab <- merge(cecum,colon[,c(1,c(7:ncol(colon)))],by="Geneid")
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
  search_table <- search_table[order(search_table$sample_name),]
  age <- search_table$age
  mouse_ID <- search_table$mouse_ID
  tissue <- search_table$tissue
  age[which(age=="3m")] <- "young"
  age[which(age=="24m")] <- "old"
  colnames(counts) <- paste0(colnames(counts),"-",age,"-",mouse_ID)
  batch <- c("batch1","batch2","batch2","batch3","batch1","batch3","batch4","batch4")
  samples <- data.frame(
    age = paste0(tissue,"_",age),
    tissue = tissue,
    batch = batch
  )
  samples$age <- factor(samples$age, levels = c("Cecum_young","Cecum_old","Colon_young","Colon_old"))
  y <- DGEList(counts=counts, group=samples$tissue)  
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  tab<-tab[keep,]
  y$samples$age <- age
  y$samples$batch <- batch
  y <- calcNormFactors(y)
  design <- model.matrix(~batch + age, data=samples)  
  colnames(design) <- gsub("tissue", "", colnames(design))  
  colnames(design) <- gsub("age", "", colnames(design))  
  y <- estimateDisp(y, design) 
  fit_tag = glmFit(y,design)
  contrast_matrix <- makeContrasts(  
    Colon_Old_vs_Young = Colon_old - Colon_young, 
    levels = design  
  )  
  lrt_colon <- glmLRT(fit_tag, contrast=contrast_matrix[,"Colon_Old_vs_Young"])  
  colon_out <-  cbind(tab[,1:6],cpm(y)[,which(tissue == "Colon") ],logCPM=lrt_colon$table$logCPM,bcv=sqrt(fit_tag$dispersion),
                      "PValue.old-young"=lrt_colon$table$PValue,"FDR.old-young"= p.adjust(lrt_colon$table$PValue,method="BH"),
                      "LogFC.old-young"=lrt_colon$table$logFC)
  colon_out$Significant <- ifelse(colon_out$`FDR.old-young` < 0.05 & abs(colon_out$`LogFC.old-young`) >= log2(1.2), 
                            ifelse(colon_out$`LogFC.old-young` > log2(1.2), "Up", "Down"), "Stable")
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  write.csv(colon_out,paste0("data/samples/colon/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  p_list[["colon"]] <- ggplot(
    colon_out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0("Colon ",antibody))+
    annotate("text", x = min(colon_out$`LogFC.old-young`), y = max(-log10(colon_out$`FDR.old-young`)), label = nrow(colon_out[which(colon_out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(colon_out$`LogFC.old-young`), y = max(-log10(colon_out$`FDR.old-young`)), label = nrow(colon_out[which(colon_out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  
  samples$age <- factor(samples$age, levels = c("Colon_young","Colon_old","Cecum_young","Cecum_old"))
  design <- model.matrix(~0+batch + age, data=samples)  
  colnames(design) <- gsub("tissue", "", colnames(design))  
  colnames(design) <- gsub("age", "", colnames(design))
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  contrast_matrix <- makeContrasts(  
    Cecum_Old_vs_Young = Cecum_old - Cecum_young, 
    levels = design  
  )  
  lrt_cecum <- glmLRT(fit_tag, contrast=contrast_matrix[,"Cecum_Old_vs_Young"]) 
  cecum_out <-  cbind(tab[,1:6],cpm(y)[,which(tissue == "Cecum") ],logCPM=lrt_cecum$table$logCPM,bcv=sqrt(fit_tag$dispersion),
                      "PValue.old-young"=lrt_cecum$table$PValue,"FDR.old-young"= p.adjust(lrt_cecum$table$PValue,method="BH"),
                      "LogFC.old-young"=lrt_cecum$table$logFC)
  cecum_out$Significant <- ifelse(cecum_out$`FDR.old-young` < 0.05 & abs(cecum_out$`LogFC.old-young`) >= log2(1.2), 
                                  ifelse(cecum_out$`LogFC.old-young` > log2(1.2), "Up", "Down"), "Stable")
  write.csv(cecum_out,paste0("data/samples/cecum/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  p_list[["cecum"]] <- ggplot(
    cecum_out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0("Cecum ",antibody))+
    annotate("text", x = min(cecum_out$`LogFC.old-young`), y = max(-log10(cecum_out$`FDR.old-young`)), label = nrow(cecum_out[which(cecum_out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(cecum_out$`LogFC.old-young`), y = max(-log10(cecum_out$`FDR.old-young`)), label = nrow(cecum_out[which(cecum_out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  return(p_list)
}
cecum_list <- list()
colon_list <- list()
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K4me3","H3K4me1","H3K27ac")
for(i in c(1:length(antibodys))){
  p_list <- diff_analysis(antibodys[i])
  colon_list[[i]] <- p_list[["colon"]]
  cecum_list[[i]] <- p_list[["cecum"]]
}
colon_combined_plot <- plot_a_list(colon_list,2,3)
cecum_combined_plot <- plot_a_list(cecum_list,2,3)
ggsave(paste0("result/colon/all_diff_volcano_plot_remove_batch_effect.png"),colon_combined_plot,width = 15,height = 10,type="cairo")
ggsave(paste0("result/cecum/all_diff_volcano_plot_remove_batch_effect.png"),cecum_combined_plot,width = 15,height = 10,type="cairo")

diff_in_peaks_volcano <- function(tissue, antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  out <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"),row.names = 1)
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    peaks <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_10kb_in_young_old_merge-W1000-G3000-E100.bed"),header = F)
  }else{
    peaks <- read.delim(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_1kb_in_young_old_merge_macs_narrowpeak.bed"),header = F)
  }
  out <- out[which(out$Geneid %in% peaks$V4),]
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(
    out, aes(x = `LogFC.old.young`, y = -log10(`FDR.old.young`))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
    annotate("text", x = min(out$`LogFC.old.young`), y = max(-log10(out$`FDR.old.young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old.young`), y = max(-log10(out$`FDR.old.young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  return(p)
}
tissue <- "colon"
p_list <- list()
for(j in c(1:length(antibodys))){
  antibody <- antibodys[j]
  p_list[[j]] <- diff_in_peaks_volcano(tissue,antibody)
}
combined_plot <- plot_a_list(p_list,2,3)
ggsave(paste0("result/",tissue,"/all_diff_volcano_plot_remove_batch_effect_in_peaks.png"),combined_plot,width = 15,height = 10,type="cairo")
