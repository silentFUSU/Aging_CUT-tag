rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
tissue <- "lung"
tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))


summary <- data.frame()
for(tissue in tissues){
  domain <- read.table(paste0(paste0("data/samples/",tissue,"/H3K27me3/peaks/edd/edd_peaks_fdr05.bed")))
  domain <- domain[,c(1:3)]
  domain <- as.data.table(domain)
  setDT(domain)
  setkey(domain,V1,V2,V3)
  
  df <- read.csv(paste0("data/samples/",tissue,"/H3K27me3/H3K27me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  df <- df[,c("Chr","Start","End","Significant")]
  df <- as.data.table(df)
  setDT(df)
  setkey(df,Chr,Start,End)
  overlaps <- foverlaps(df, domain, type = "any", nomatch = 0L)  
  t_summary <- as.data.frame(table(overlaps$Significant))
  t_summary$percentage <- t_summary$Freq/sum(t_summary$Freq)*100
  t_summary$tissue <- tissue_label_change(tissue)
  t_summary <- t_summary[,c("Var1","percentage","tissue")]
  summary <- rbind(summary,t_summary)
}

summary_sort <- summary[which(summary$Var1=="Up"),]
summary_sort <- summary_sort[order(summary_sort$percentage),]
tissue_sort <- summary_sort$tissue
summary$Var1 <- factor(summary$Var1,levels=c("Up","Stable","Down"))
summary$tissue <- factor(summary$tissue, levels=tissue_sort)
color <- setNames(c("#009980","#E69900","#838B8B"),c("Up","Stable","Down"))
ggplot(summary, aes(x = tissue, y = percentage, fill = Var1)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  # geom_text(data = to_plot,
  #           aes(label = Freq, y = position),
  #           color = "black", size = 5, vjust = 0.5) +
  ylab("Proportion")+
  ggtitle("H3K27me3 bin change in age domains")


summary <- data.frame()
for(tissue in tissues){
  domain <- read.table(paste0(paste0("data/samples/",tissue,"/H3K27me3/peaks/edd/edd_peaks_fdr05.bed")))
  domain <- domain[,c(1:3)]
  domain <- as.data.table(domain)
  setDT(domain)
  setkey(domain,V1,V2,V3)
  
  df <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  df <- df[,c("Chr","Start","End","Significant")]
  df <- as.data.table(df)
  setDT(df)
  setkey(df,Chr,Start,End)
  overlaps <- foverlaps(df, domain, type = "any", nomatch = 0L)  
  t_summary <- as.data.frame(table(overlaps$Significant))
  t_summary$percentage <- t_summary$Freq/sum(t_summary$Freq)*100
  t_summary$tissue <- tissue_label_change(tissue)
  t_summary <- t_summary[,c("Var1","percentage","tissue")]
  summary <- rbind(summary,t_summary)
}

# summary_sort <- summary[which(summary$Var1=="Down"),]
# summary_sort <- summary_sort[order(summary_sort$percentage),]
# tissue_sort <- summary_sort$tissue
summary$Var1 <- factor(summary$Var1,levels=c("Up","Stable","Down"))
summary$tissue <- factor(summary$tissue, levels=tissue_sort)
color <- setNames(c("#009980","#E69900","#838B8B"),c("Up","Stable","Down"))
ggplot(summary, aes(x = tissue, y = percentage, fill = Var1)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  # geom_text(data = to_plot,
  #           aes(label = Freq, y = position),
  #           color = "black", size = 5, vjust = 0.5) +
  ylab("Proportion")+
  ggtitle("H3K9me3 bin change in H3K27me3 age domains")



antibody <- "H3K27me3"
diff_bin_level_remove_batch_effect <- function(tissue,antibody){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  tab =  read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_edd_peaks_fdr05.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  counts <- counts[,search_table$sample_name]
  search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
  search_table$age <- factor(search_table$age, levels = c("3m","24m"))
  
  age <- as.character(search_table$age)
  batch <- as.character(search_table$batch)
  mouse_ID <- search_table$mouse_ID
  age[which(age=="3m")] <- "young"
  age[which(age=="24m")] <- "old"
  colnames(counts) <- paste0(colnames(counts),"-",age,"-",mouse_ID,"-",batch)
  
  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$year <- age
  y$samples$year <- factor(y$samples$year,c("young","old"))
  y$samples$batch <- search_table$batch
  
  y <- calcNormFactors(y)
  design <- model.matrix(~batch+year, y$samples)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = which(colnames(design) == "yearold"))
  tab<-tab[keep,]
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
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
    ggtitle(paste0(tissue_label_change(tissue)," H3K27me3 age-domains"))+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  write.csv(out,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_edd_peaks_fdr05_diff_after_remove_batch_effect.csv"),row.names = F)
}
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
antibody <- "H3K27me3"
p_list <- list()
for(tissue in tissues){
  p_list[[tissue]] <-   diff_bin_level_remove_batch_effect(tissue,antibody)
}