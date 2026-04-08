.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(corrplot)
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
antibodys <- c("H4K20me3","H3K9me2","H2AK119ub")
for(antibody in antibodys){
  tab <- read.delim(paste0("data/samples/lung/",antibody,"/",antibody,"_10kb_bins.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  age <- c("young","old")
  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  cpm <- cpm(y)
  colnames(cpm) <- paste(rep(antibody,2),colnames(cpm),age,sep = "-")
  p <- ggplot(cpm, aes(x = log10(cpm[,1]), y = log10(cpm[,2]))) +  
    geom_point(color="grey") +  
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") + 
    labs(  
      title = antibody,  
    ) +                 
    xlab(paste(colnames(cpm)[1])) +
    ylab(paste(colnames(cpm)[2])) +
    theme_minimal()   
  # print(p)
  
  y$samples$year <- age
  y$samples$year <- factor(y$samples$year,c("young","old"))
  y <- calcNormFactors(y)
  design <- model.matrix(~year, y$samples)
  bcv <- 0.1
  fit_tag = glmFit(y,design,dispersion = bcv^2)
  lrt = glmLRT(fit_tag, coef = 2)
  t_tab<-tab[keep,]
  out = cbind(t_tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
              "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC.old-young"=lrt$table$logFC)
  out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= log2(1.2), 
                            ifelse(out$`LogFC.old-young` > log2(1.2), "Up", "Down"), "Stable")
  
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  
  p <- ggplot(
    out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0("Lung ",antibody))+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  print(p)
  }

tab <- read.delim(paste0("data/raw_data/20250120_LLX_CUTTag/10kb_bins.counts"),skip=1)
counts = tab[,c(7:ncol(tab))]
rownames(counts)= tab$Geneid
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
age <- factor(c("young","old"),c("young","old"))
antibodys <- c("H4K20me3","H3K79me3","H3K9me2","H2A_K119ub")
for(i in c(1:4)){
  t_counts <- counts[,c(i,i+4)]  
  y= DGEList(counts=t_counts,group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$year <- age
  y$samples$year <- factor(y$samples$year,c("young","old"))
  y <- calcNormFactors(y)
  design <- model.matrix(~year, y$samples)
  bcv <- 0.1
  fit_tag = glmFit(y,design,dispersion = bcv^2)
  lrt = glmLRT(fit_tag, coef = 2)
  t_tab<-tab[keep,]
  out = cbind(t_tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
              "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC.old-young"=lrt$table$logFC)
  out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= log2(1.2), 
                            ifelse(out$`LogFC.old-young` > log2(1.2), "Up", "Down"), "Stable")
  
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))

  p <- ggplot(
    out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0("Lung ",antibodys[i]))+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  print(p)
}


antibodys <- c("H3K27me3","H3K9me3","H4K20me3","H3K9me2","H2AK119ub")
df_summary <- data.frame()
# df <- read.csv("data/samples/lung/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv")
# df <- read.table("data/samples/lung/H3K27me3_H3K9me3_intersect/10kb_all_significant_second_quadrant_after_remove_batch_effect.bed")
# regions <- df$V4
for(antibody in antibodys){
  if(antibody %in% c("H3K27me3","H3K9me3")){
    df <- read.csv(paste0("data/samples/lung/",antibody,"/",antibody,"_10kb_bins_diff_after_remove_batch_effect.csv"))
    # df <- df[which(df$Geneid %in% regions),c("Geneid","LogFC.old.young")]
    df <- df[,c("Geneid","LogFC.old.young")]
    colnames(df)[which(colnames(df) == "LogFC.old.young")] <- antibody
    if(nrow(df_summary)==0){
      df_summary <- df
    }else{
      df_summary <- merge(df_summary,df,by="Geneid")
    }
  }else{
    tab <- read.delim(paste0("data/samples/lung/",antibody,"/",antibody,"_10kb_bins.counts"),skip=1)
    counts = tab[,c(7:ncol(tab))]
    rownames(counts)= tab$Geneid
    pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
    colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
    counts$LogFC.old.young <- log2(counts[,2]/counts[,1])
    counts$Geneid <- rownames(counts)
    df <- counts[,c("Geneid","LogFC.old.young")]
    colnames(df)[which(colnames(df) == "LogFC.old.young")] <- antibody
    if(nrow(df_summary)==0){
      df_summary <- df
    }else{
      df_summary <- merge(df_summary,df,by="Geneid")
    }
  }
}
rownames(df_summary) <- df_summary$Geneid
df_summary <- df_summary[,-1]
df_summary <- df_summary[order(df_summary$H3K27me3),]
df_summary <- df_summary[!apply(df_summary, 1, function(row) any(row == Inf | row == -Inf)), ]
pheatmap::pheatmap(df_summary,show_rownames = F,cluster_rows = F,breaks = seq(-2, 2, length.out = 101))

cor <- cor(df_summary)
pheatmap::pheatmap(cor,breaks = seq(-1, 1, length.out = 101))

df <- df[which(df$Geneid %in% out$Geneid[which(out$Significant=="Up")]),]
