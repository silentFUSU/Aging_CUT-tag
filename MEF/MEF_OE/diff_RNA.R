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
# conditions <- c("MEF_Bmi1", "MEF_Cbx2", "MEF_Cbx7")
conditions <- c("MEF_mEzh2", "MEF_hEzh2", "MEF_mCbx8")
p_list <- list()

for(condition in conditions){
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  if(condition %in% c("MEF_mEzh2", "MEF_hEzh2")){
    tab = read.delim(paste0("data/samples/RNA/MEF_OE_RNA/MEF_Vector2_",condition,"_merge.counts"),skip=1)
  }else{
    tab = read.delim(paste0("data/samples/RNA/MEF_OE_RNA/MEF_Vector_",condition,"_merge.counts"),skip=1)
  }

  tab <- tab[!grepl("chrY", tab$Chr), ]
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  new_names <- sub("^.*\\.bam\\.([A-Za-z]+\\d+\\.\\d+)(?:[_.].*)?$", "\\1", colnames(counts))
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
  # keep = which(rowSums(cpm(y)>0)>=2)
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
  out$Significant <- ifelse(out$`FDR.oe-vec` < 0.05 & abs(out$`LogFC.oe-vec`) >= 0, 
                            ifelse(out$`LogFC.oe-vec` > 0, "Up", "Down"), "Stable")
  if(condition %in% c("MEF_mEzh2", "MEF_hEzh2")){
    write.csv(out,paste0("data/samples/RNA/MEF_OE_RNA/MEF_Vector2_",condition,"_diff_expression_gene_strict_filter_bar.csv"),row.names = F)
  }else{
    write.csv(out,paste0("data/samples/RNA/MEF_OE_RNA/MEF_Vector_",condition,"_diff_expression_gene_strict_filter_bar.csv"),row.names = F)
  }
  
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  
  top_genes <- out[which(out$Significant=="Up"),] %>%  
    arrange(desc(`LogFC.oe-vec`)) %>%  
    head(10)  
  bottom_genes <- out[which(out$Significant=="Down"),] %>%  
    arrange(`LogFC.oe-vec`) %>%  
    head(10)  
  highlight_genes <- rbind(top_genes, bottom_genes) 
  highlight_genes$Gene <- rownames(highlight_genes)
  p_list[[paste0(condition)]] <- ggplot(
    out, aes(x = `LogFC.oe-vec`, y = -log10(`FDR.oe-vec`))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(condition))+
    annotate("text", x = min(out$`LogFC.oe-vec`), y = max(-log10(out$`FDR.oe-vec`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.oe-vec`), y = max(-log10(out$`FDR.oe-vec`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)+
    geom_text_repel(data = highlight_genes, aes(`LogFC.oe-vec`,-log10(`FDR.oe-vec`), label = Geneid), max.overlaps=100,
                   size = 5, 
                   nudge_y = 0.2)
}

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_p <- plot_a_list(p_list,no_of_rows=1,no_of_cols=3)
# combined_p
ggsave("result/MEF_OE/all_diff_gene_strict_filter_bar2.png",combined_p,width=15,height=5)

