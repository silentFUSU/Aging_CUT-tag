rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(stringr)
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
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}

tissue <- "lung"
peak_preprocess_peak_level_remove_batch_effect <- function(tissue,antibody){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_recursion.counts"),skip=1)
  if(tissue %in% c("mammarygland","ovary","uterus")){
    tab <- tab[-which(tab$Chr=="chrY"),]
  }
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|DYQ[0-9]+).*"
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

  write.csv(out,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_young_old_merge-W",window_size,"-G",gap_size,"-E100_recursion_diff_after_remove_batch_effect.csv"),row.names = F)
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(
      out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
      geom_point(aes(color = Significant, size = Length),alpha=0.2) +
      scale_color_manual(values = colour) +
      scale_size_continuous(range = c(1, 5)) +
      geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
      labs(x="log2(fold change)",
           y="-log10 (FDR)") +
      theme_bw()+
      theme(text = element_text(size = 20))+
      ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
      annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
      annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  out <- out[which(out$Length > 200000),]
  p <- ggplot(
    out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
    geom_point(aes(color = Significant, size = Length),alpha=0.2) +
    scale_color_manual(values = colour) +
    scale_size_continuous(range = c(1, 5)) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (FDR)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  return(p)
}

p_list <- list()
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))
window_size="5000"
gap_size="10000"
antibody <- "H3K9me3"
for(j in c(1:length(tissues))){
  tissue <- tissues[j]
  p_list[[j]] <- peak_preprocess_peak_level_remove_batch_effect(tissue,antibody)
}
combined_plot <- plot_a_list(p_list,4,7)
ggsave(paste0("result/all/diff/",antibody,"/all_tissues_diff_-W",window_size,"-G",gap_size,"-E100_recursion_peaks_volcano_plot_remove_batch_effect.png"),combined_plot,width = 40,height = 18,type="cairo")



