rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(stringr)
tissue <- "ileum"
antibody <- "H3K9me3"
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
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
    }
  }
  return(tissue_label)
}
label_shuffling_diff_analysis <- function(tissue,antibody){
  if(antibody %in% c("H3K9me3","H3K27me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  if (tissue=="ovary"){
    colnames(counts) = c("old_1","young_1","old_2","young_2")
    group =c("label2","label1","label1","label2")
  }else{
    colnames(counts) = c("young_1","old_1","young_2","old_2")
    group =c("label1","label2","label2","label1")
  }
  
  y= DGEList(counts=counts,group=group)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$batch <- rep(c(rep("batch1", 2), rep("batch2", 2)), 1)
  y <- calcNormFactors(y)
  design <- model.matrix(~batch+group, y$samples)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = 3)
  tab<-tab[keep,]
  
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
              "PValue"=lrt$table$PValue,"FDR"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC"=lrt$table$logFC)
  
  out$Significant <- ifelse(out$FDR < 0.05 & abs(out$LogFC) >= 0, 
                            ifelse(out$LogFC > 0, "Up", "Down"), "Stable")
  # write.csv(out,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_label_shuffling_diff.csv"),row.names = F)
  color <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(
    out, aes(x = LogFC, y = -log10(FDR))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = color) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (FDR)") +
    theme_bw()+
    ggtitle(tissue_label_change(tissue))+
    theme(text = element_text(size = 20))+
    annotate("text", x = min(out$LogFC), y = max(-log10(out$FDR)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$LogFC), y = max(-log10(out$FDR)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  print(p)
  return(p)
}
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder","CB","jejunum","uterus","ovary")
antibodys <- c("H3K27me3","H3K36me3","H3K9me3","H3K27ac","H3K4me3","H3K4me1")
for(antibody in antibodys){
  p_list <- list()
  for(i in c(1:length(tissues))){
    tissue <- tissues[i]
    p_list[[i]] <- label_shuffling_diff_analysis(tissue,antibody)
  }
  combined_plot <- plot_a_list(p_list,4,6)
  ggsave(paste0("result/all/diff/",antibody,"/label_shuffling_diff_analysis.png"),combined_plot,width = 40,height = 25,type="cairo")
}
