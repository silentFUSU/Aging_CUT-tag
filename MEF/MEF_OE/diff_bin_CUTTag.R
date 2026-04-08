rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
antibodys <- c("H3K27me3","H2AK119ub1")
conditions <- c("MEF_Bmi1", "MEF_Cbx2", "MEF_Cbx7")
p_list <- list()
for(antibody in antibodys){
  for(condition in conditions){
    search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
    tab = read.delim(paste0("data/samples/MEF_OE/",antibody,"/",antibody,"_MEF_Vector_",condition,"_10kb.counts"),skip=1)
    tab <- tab[!grepl("chrY", tab$Chr), ]
    counts = tab[,c(7:ncol(tab))]
    rownames(counts)= tab$Geneid
    new_names <- sub("^.*\\.bam\\.([A-Za-z]+\\d+\\.\\d+)\\.bam.*$", "\\1", colnames(counts))
    new_names <- gsub("\\.", "-", new_names)
    colnames(counts) <- new_names
    search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
    counts <- counts[,search_table$sample_name]
    search_table$age <- factor(search_table$age, levels = c("vec","oe"))
    age <- as.character(search_table$age)
    batch <- as.character(search_table$batch)
    mouse_ID <- search_table$mouse_ID
    
    colnames(counts) <- paste0(colnames(counts),"-",age,"-",mouse_ID,"-",batch)
    y= DGEList(counts=counts,group=age)
    keep = which(rowSums(cpm(y)>1)>=2)
    y = y[keep,]
    y$samples$year <- age
    y$samples$year <- factor(y$samples$year,c("vec","oe"))
    y$samples$batch <- search_table$batch
    y <- calcNormFactors(y)
    if(length(unique(y$samples$batch))==1){
      design <- model.matrix(~year, y$samples)
    }else{
      design <- model.matrix(~batch+year, y$samples)
    }
    y<-estimateCommonDisp(y)
    y<-estimateGLMTagwiseDisp(y,design)
    fit_tag = glmFit(y,design)
    lrt = glmLRT(fit_tag, coef = which(colnames(design) == "yearoe"))
    tab<-tab[keep,]
    
    out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
                "PValue.oe-vec"=lrt$table$PValue,"FDR.oe-vec"= p.adjust(lrt$table$PValue,method="BH"),
                "LogFC.oe-vec"=lrt$table$logFC)
    out$Significant <- ifelse(out$`FDR.oe-vec` < 0.05 & abs(out$`LogFC.oe-vec`) >= log2(1.2), 
                              ifelse(out$`LogFC.oe-vec` > log2(1.2), "Up", "Down"), "Stable")
    write.csv(out,paste0("data/samples/MEF_OE/",antibody,"/",condition,"/",antibody,"_",condition,"_10kb_bins_diff_after_remove_batch_effect.csv"),row.names = F)
    colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
    p_list[[paste0(condition,"-",antibody)]] <- ggplot(
      out, aes(x = `LogFC.oe-vec`, y = -log10(`FDR.oe-vec`))) +
      geom_point(aes(color = Significant), size=2) +
      scale_color_manual(values = colour) +
      geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
      labs(x="log2(foe change)",
           y="-log10 (p-value)") +
      theme_bw()+
      theme(text = element_text(size = 20),legend.position = "none")+
      ggtitle(paste0(condition," ",antibody))+
      annotate("text", x = min(out$`LogFC.oe-vec`), y = max(-log10(out$`FDR.oe-vec`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
      annotate("text", x = max(out$`LogFC.oe-vec`), y = max(-log10(out$`FDR.oe-vec`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
    }
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_p <- plot_a_list(p_list,no_of_rows=2,no_of_cols=3)
combined_p
ggsave("result/MEF_OE/all_diff_bin.png",combined_p,width=15,height=10)
