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
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
Gene_TSS_diff_remove_bath_effect <- function(tissue,antibody){
  p_list <- list()
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_gene_TSS_100kb.counts"),skip=1)
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
  
  write.csv(out,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_gene_TSS_100kb_diff_after_remove_batch_effect.csv"),row.names = F)
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p_list[[1]] <- ggplot(
    out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W1000-G3000-E100.bed"))
  peaks <- as.data.table(peaks)
  setDT(peaks)
  setkey(peaks,V1,V2,V3)

  gene_tss <- out[,c(1:4)]
  gene_tss <- as.data.table(gene_tss)
  gene_tss$Start <- as.numeric(gene_tss$Start)
  setDT(gene_tss)
  setkey(gene_tss,Chr,Start,End)
  overlaps <- foverlaps(gene_tss,peaks, type = "any", nomatch = 0L)
  out_inpeak <- out[which(out$Geneid %in% overlaps$Geneid),]
  p_list[[2]] <- ggplot(
    out_inpeak, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))+
    annotate("text", x = min(out_inpeak$`LogFC.old-young`), y = max(-log10(out_inpeak$`FDR.old-young`)), label = nrow(out_inpeak[which(out_inpeak$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out_inpeak$`LogFC.old-young`), y = max(-log10(out_inpeak$`FDR.old-young`)), label = nrow(out_inpeak[which(out_inpeak$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  out <- out[which(out$Significant != "Stable"),]
  rna <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_filter_bar.csv"))
  # rna <- rna[which(rna$Significant!="Stable"),]
  colnames(rna)[1] <- "Geneid"
  to_plot <- merge(rna[,c("Geneid","logFC","Significant")],out[,c("Geneid","LogFC.old-young")],by="Geneid")
  colnames(to_plot) <- c("Geneid","RNA_logFC","RNA_Significant","histone_logFC")
  
  top_genes <- to_plot[which(to_plot$RNA_Significant=="Up"),] %>%  
    arrange(desc(RNA_logFC)) %>%  
    head(10)  
  bottom_genes <- to_plot[which(to_plot$RNA_Significant=="Down"),] %>%  
    arrange(RNA_logFC) %>%  
    head(10)  
  highlight_genes <- rbind(top_genes, bottom_genes)  
  x_range <- range(to_plot$RNA_logFC, na.rm = TRUE)  
  y_range <- range(to_plot$histone_logFC, na.rm = TRUE)  
  x_pos_right <- x_range[2] * 0.9    
  x_pos_left <- x_range[1] * 0.9   
  y_pos_top <- y_range[2] * 0.9    
  y_pos_bottom <- y_range[1] * 0.9 
  
  p_list[[3]] <- ggplot() +
    geom_point(data=to_plot, mapping=aes(RNA_logFC, histone_logFC),color = "grey",alpha=0.5)+
    geom_point(data=to_plot[which(to_plot$RNA_Significant=="Up"),], mapping=aes(RNA_logFC, histone_logFC),color = "red") +
    geom_point(data=to_plot[which(to_plot$RNA_Significant=="Down"),], mapping=aes(RNA_logFC, histone_logFC),color = "blue") +
    labs(x="RNA log2(Fold Change)",
         y=paste0(antibody," log2(Fold Change)")) +
    geom_point(color="grey")+
    theme(legend.position = "bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme_bw()+theme(text = element_text(size = 18))+
    # ggtitle(tissue_label_change(tissue),paste0("Chi-squared test ",chi_label_summary[which(chi_label_summary$tissue==tissue_label_change(tissue)),"label"]))+
    ggtitle(tissue_label_change(tissue))+
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_text_repel(data = highlight_genes, aes(RNA_logFC, histone_logFC, label = Geneid),   
                    size = 5, # 字体大小  
                    nudge_y = 0.2)+
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$RNA_logFC > 0 & to_plot$histone_logFC > 0 & to_plot$RNA_Significant == "Up"), ])),  
             x = x_pos_right, y = y_pos_top, colour = "#00b8a9", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$RNA_logFC < 0 & to_plot$histone_logFC < 0 & to_plot$RNA_Significant == "Down"), ])),  
             x = x_pos_left, y = y_pos_bottom, colour = "#ff9a00", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$RNA_logFC < 0 & to_plot$histone_logFC > 0 & to_plot$RNA_Significant == "Down"), ])),  
             x = x_pos_left, y = y_pos_top, colour = "#f6416c", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$RNA_logFC > 0 & to_plot$histone_logFC < 0 & to_plot$RNA_Significant == "Up"), ])),  
             x = x_pos_right, y = y_pos_bottom, colour = "#48466d", size = 5) 
  return(p_list)
  }
tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
p_list <- list()
volcano_plot <- list()
volcano_plot_in_peaks <- list()
dot_plot <- list()
for(tissue in tissues){
  p_list <- Gene_TSS_diff_remove_bath_effect(tissue, "H3K9me3")
  volcano_plot[[tissue]] <- p_list[[1]]
  volcano_plot_in_peaks[[tissue]] <- p_list[[2]]
  dot_plot[[tissue]] <- p_list[[3]]
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
volcano_plot_combined <- plot_a_list(volcano_plot,no_of_rows = 4,no_of_cols = 7)
ggsave(paste0("result/RNA/histone_relationship_with_RNA/",antibody,"/",antibody,"_relationship_with_RNA_gene_TSS_diff_volcano_plot.png"),volcano_plot_combined,width = 42,height = 24,type="cairo")
volcano_plot_in_peaks_combined <- plot_a_list(volcano_plot_in_peaks,no_of_rows = 4,no_of_cols = 7)
ggsave(paste0("result/RNA/histone_relationship_with_RNA/",antibody,"/",antibody,"_relationship_with_RNA_gene_TSS_diff_overlap_with_peak_volcano_plot.png"),volcano_plot_in_peaks_combined,width = 42,height = 24,type="cairo")
dot_plot_combined <- plot_a_list(dot_plot,no_of_rows = 4,no_of_cols = 7)
ggsave(paste0("result/RNA/histone_relationship_with_RNA/",antibody,"/",antibody,"_relationship_with_RNA_gene_TSS_diff_relationship_with_gene_dot_plot.png"),dot_plot_combined,width = 42,height = 24,type="cairo")

