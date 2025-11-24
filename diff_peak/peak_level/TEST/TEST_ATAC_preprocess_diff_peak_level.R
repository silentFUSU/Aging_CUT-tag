rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
bin_size <-"1kb"
antibody <- "ATAC"
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
peak_preprocess <- function(tissue){
  tab = read.delim(paste0("data/samples/ATAC/",tissue,"/ATAC/ATAC_macs_young_old_narrowpeak.counts"),skip=1)
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(tab)[7:length(tab)] <-  gsub(pattern, "\\1", colnames(tab)[7:length(tab)])
  counts <- tab[7:length(tab)]
  search_table <- read.csv("data/samples/all/ATAC_search_table.csv")
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
  search_table <- search_table[order(search_table$sample_name),]
  if(tissue == "testis"){
    counts <- counts[,search_table$sample_name]
  }
  age <- search_table$age
  mouse_ID <- search_table$mouse_ID
  age[which(age=="3m")] <- "young"
  age[which(age=="24m")] <- "old"
  colnames(counts) <- paste0(colnames(counts),"-",age,"-",mouse_ID)

  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$group <- factor(y$samples$group,c("young","old"))
  y <- calcNormFactors(y)
  design <- model.matrix(~group, y$samples)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = 2)
  tab<-tab[keep,]
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
              "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC.old-young"=lrt$table$logFC)
  out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                            ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
  write.csv(out,paste0("data/samples/ATAC/",tissue,"/ATAC_peaks_diff.csv"),row.names = F)
  outup <- out[which(out$Significant=="Up"),]
  outdown <- out[which(out$Significant=="Down"),]
  write.table(outdown[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/ATAC/",tissue,"/ATAC/bed/ATAC_peaks_diff_down.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(outup[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/ATAC/",tissue,"/ATAC/bed/ATAC_peaks_diff_up.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(
    out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (fdr)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," ATAC"))+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  return(p)
}
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
p_list <- list()
for(i in c(1:length(tissues))){
  p_list[[i]] <- peak_preprocess(tissues[i])
}
combined_plot <- plot_a_list(p_list,no_of_rows = 4,no_of_cols = 7)
ggsave("result/all/ATAC/all_tissues_peaks_diff.png",combined_plot,width = 35,height = 20,type="cairo")


split_diff_peaks_bed <- function(tissue){
  df <- read.csv(paste0("data/samples/ATAC/",tissue,"/ATAC_peaks_diff.csv"))
  outup <- df[which(df$Significant=="Up"),]
  outdown <- df[which(df$Significant=="Down"),]
  outstable <- df[which(df$Significant=="Stable"),]
  write.table(outdown[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/ATAC/",tissue,"/ATAC/bed/ATAC_peaks_diff_down.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(outup[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/ATAC/",tissue,"/ATAC/bed/ATAC_peaks_diff_up.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(outstable[,c("Chr","Start","End","Geneid")], file=paste0("data/samples/ATAC/",tissue,"/ATAC/bed/ATAC_peaks_diff_stable.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
for(tissue in tissues){
  split_diff_peaks_bed(tissue)
}
