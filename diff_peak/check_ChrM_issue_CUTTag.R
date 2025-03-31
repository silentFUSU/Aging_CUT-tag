rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
tissues <- c("uterus","mammarygland")
antibodys <- c("H3K9me3","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K27ac")
issue_mouse <- c("100","105","125","135","136","138","139","140","203","205","215","235","212","213","224","225","226","230","233","245","246","250")
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
chrM_issue_diff <- function(tissue,antibody){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
  search_table <- search_table[which(search_table$age=="3m"),]
  counts <- counts[,search_table$sample_name]
  mouse <- as.character(search_table$mouse_ID)
  
  y= DGEList(counts=counts,group=mouse)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$group <- factor(y$samples$group, levels=c(mouse[-which(mouse %in% wrong_mouse)],mouse[which(mouse %in% wrong_mouse)]))
  design <- model.matrix(~group, y$samples)
  y <- calcNormFactors(y)
  bcv <- 0.2
  fit_tag = glmFit(y,design,dispersion = bcv^2)
  lrt = glmLRT(fit_tag, coef = 2)
  t_tab<-tab[keep,]
  out = cbind(t_tab[,c(1:6)],cpm(y),lrt$table, "fdr"=p.adjust(lrt$table$PValue,method="BH"))
  out$Significant <- ifelse(out$fdr< 0.05 & abs(out$logFC) >= 0, 
                            ifelse(out$logFC > 0, "Up", "Down"), "Stable")
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  
  p <- ggplot(
    out, aes(x = logFC, y = -log10(fdr))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody," ",mouse[which(mouse %in% wrong_mouse)]," vs ",mouse[-which(mouse %in% wrong_mouse)]))+
    annotate("text", x = min(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  print(p)
  dir.create("data/chrM_inssue/CUTTag/")
  write.csv(out,paste0("data/chrM_inssue/CUTTag/",tissue,"_",antibody,"_",mouse[which(mouse %in% wrong_mouse)],"_vs_",mouse[-which(mouse %in% wrong_mouse)],"_bcv02.csv"))
  return(p)
  }
p_list <- list()
for(antibody in antibodys){
  p_list[[antibody]] <- chrM_issue_diff(tissue,antibody)
}
combined_plot <- plot_a_list(p_list,2,3)
ggsave(paste0("data/chrM_inssue/CUTTag/",tissue,"_",mouse[which(mouse %in% wrong_mouse)],"_vs_",mouse[-which(mouse %in% wrong_mouse)],"_bcv02.png"),combined_plot,width = 15,height = 10,type="cairo")

chrM_issue_diff_with_rep <- function(tissue,antibody){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
  search_table <- search_table[which(search_table$age=="3m"),]
  counts <- counts[,search_table$sample_name]
  mouse <- as.character(search_table$mouse_ID)
  mouse[-which(mouse %in% wrong_mouse)] <- "other"
  y= DGEList(counts=counts,group=mouse)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$group <- factor(y$samples$group, levels=c("other",mouse[which(mouse %in% wrong_mouse)]))
  design <- model.matrix(~group, y$samples)
  y <- calcNormFactors(y)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = ncol(design))
  tab<-tab[keep,]
  out = cbind(tab[,c(1:6)],cpm(y),lrt$table, "fdr"=p.adjust(lrt$table$PValue,method="BH"))
  out$Significant <- ifelse(out$fdr< 0.05 & abs(out$logFC) >= 0, 
                            ifelse(out$logFC > 0, "Up", "Down"), "Stable")
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  
  p <- ggplot(
    out, aes(x = logFC, y = -log10(fdr))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody," ",mouse[which(mouse %in% wrong_mouse)]," vs ",mouse[-which(mouse %in% wrong_mouse)]))+
    annotate("text", x = min(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  print(p)
  dir.create("data/chrM_inssue/CUTTag/")
  write.csv(out,paste0("data/chrM_inssue/CUTTag/",tissue,"_",antibody,"_",mouse[which(mouse %in% wrong_mouse)],"_vs_",mouse[-which(mouse %in% wrong_mouse)],"_bcv02.csv"))
  return(p)
}

tissue <- "muscle"
antibody <- "H3K36me3"
chrM_issue_diff_with_rep_thymus_muscle <- function(tissue,antibody){
  search_table <- read.csv(paste0("data/raw_data/",tissue,"_cut_tag_test/H3K36me3/search_table.csv"))
  tab = read.delim(paste0("data/raw_data/",tissue,"_cut_tag_test/H3K36me3/10kb_bins.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  search_table <- search_table[which(search_table$tissue==tissue & search_table$antibody==antibody),]
  search_table$chrM <- "normal"
  search_table$chrM[which(search_table$mouse_ID %in% wrong_mouse)] <- "mut"
  search_table <- search_table[which(search_table$condition=="paired"),]
  counts <- counts[,search_table$sample_name]
  y= DGEList(counts=counts)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  logCPMs <- cpm(y, log = TRUE)
  pca <- prcomp(t(logCPMs))
  to_plot <- data.frame(pca$x)
  to_plot$rownames <- rownames(to_plot)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  table <- search_table[which(search_table$sample_name  %in% to_plot$rownames),c(3:5,7)]
  to_plot$rownames <- rownames(to_plot)
  to_plot <- merge(to_plot,table,by.x="rownames",by.y = "sample_name")
  to_plot$label <- paste0(to_plot$rownames,"-",to_plot$mouse_ID,"-",to_plot$age,"-",to_plot$chrM)
  to_plot$age <- factor(to_plot$age,levels=c("3m","24m"))
  ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
    geom_text_repel(  
      data = to_plot,  
      aes(x = PC1, y = PC2, label = label, color = age),  
      size = 5,  
      box.padding = unit(0.35, "lines"),  
      point.padding = unit(0.3, "lines")  
    ) +
    ggtitle(paste(tissue_label_change(tissue), antibody))
  
  search_table <- search_table[which(search_table$age=="24m"),]
  
  
  
  counts <- counts[,search_table$sample_name]
  y= DGEList(counts=counts,group=search_table$chrM)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$group <- factor(y$samples$group,levels=c("normal","mut"))
  design <- model.matrix(~group, y$samples)

  y <- calcNormFactors(y)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = ncol(design))
  tab<-tab[keep,]
  out = cbind(tab[,c(1:6)],cpm(y),lrt$table, "fdr"=p.adjust(lrt$table$PValue,method="BH"))
  out$Significant <- ifelse(out$fdr< 0.05 & abs(out$logFC) >= 0, 
                            ifelse(out$logFC > 0, "Up", "Down"), "Stable")
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(
    out, aes(x = logFC, y = -log10(fdr))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 13),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody," chrM mutant vs normal"))+
    annotate("text", x = min(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  print(p)
  }
