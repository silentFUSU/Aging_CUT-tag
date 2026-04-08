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

antibodys <- c("H3K27me3","H2AK119ub1")
conditions <- c("MEF_Bmi1", "MEF_Cbx2", "MEF_Cbx7")
p_list <- list()
second_p_list <- list()
fourth_p_list <- list()
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
for(antibody in antibodys){
  for(condition in conditions){
    search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
    tab = read.delim(paste0("data/samples/MEF_OE/",antibody,"/",antibody,"_MEF_Vector_",condition,"_gene_TSS_10kb.counts"),skip=1)
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
    write.csv(out,paste0("data/samples/MEF_OE/",antibody,"/",condition,"/",antibody,"_",condition,"_gene_TSS_10kb_diff_after_remove_batch_effect.csv"),row.names = F)
    peaks <- read.table(paste0("data/samples/MEF_OE/",antibody,"/",condition,"/bed/",antibody,"_",condition,"_merge-W5000-G10000-E100.bed"))
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
    
    rna <- read.csv(paste0("data/samples/RNA/MEF_OE_RNA/MEF_Vector_",condition,"_diff_expression_gene.csv"))
    out_inpeak <- out_inpeak[which(out_inpeak$Significant != "Stable"),]
    to_plot <- merge(rna[,c("Geneid","LogFC.oe.vec","Significant")],out_inpeak[,c("Geneid","LogFC.oe-vec")],by="Geneid")
    colnames(to_plot) <- c("Geneid","RNA_logFC","RNA_Significant","histone_logFC")
    
    top_genes <- to_plot[which(to_plot$RNA_Significant=="Up"),] %>%  
      arrange(desc(RNA_logFC)) %>%  
      head(5)  
    bottom_genes <- to_plot[which(to_plot$RNA_Significant=="Down"),] %>%  
      arrange(RNA_logFC) %>%  
      head(5)  
    highlight_genes <- rbind(top_genes, bottom_genes)  
    highlight_genes <- to_plot[which(to_plot$Geneid=="Cdkn2a"),]
    x_range <- range(to_plot$RNA_logFC, na.rm = TRUE)  
    y_range <- range(to_plot$histone_logFC, na.rm = TRUE)  
    x_pos_right <- x_range[2] * 0.9    
    x_pos_left <- x_range[1] * 0.9   
    y_pos_top <- y_range[2] * 0.9    
    y_pos_bottom <- y_range[1] * 0.9 
    p_list[[paste0(condition,"_",antibody)]] <- ggplot() +
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
      ggtitle(paste0(condition,"_",antibody))+
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
    second <- to_plot$Geneid[which(to_plot$RNA_logFC < 0 & to_plot$histone_logFC > 0 & to_plot$RNA_Significant == "Down")]
    second <- bitr(second,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
    second_GO <- enrichGO( second$ENTREZID,#GO富集分析
                                OrgDb = GO_database,
                                keyType = "ENTREZID",#设定读取的gene ID类型
                                ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                                pvalueCutoff = 0.05,#设定p值阈值
                                qvalueCutoff = 0.05,#设定q值阈值
                                readable = T)
    write.csv(second_GO@result,"result/MEF_OE/Cbx7_OE_H3K27me3_RNA_second.csv")
    second_p_list[[paste0(antibody,"_",condition,"_second")]] <- barplot(second_GO,label_format = 50,showCategory = 10)+ggtitle(paste0(condition," second"))
    fourth <- to_plot$Geneid[which(to_plot$RNA_logFC > 0 & to_plot$histone_logFC < 0 & to_plot$RNA_Significant == "Up")]
    fourth <- bitr(fourth,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
    fourth_GO <- enrichGO( fourth$ENTREZID,#GO富集分析
                           OrgDb = GO_database,
                           keyType = "ENTREZID",#设定读取的gene ID类型
                           ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                           pvalueCutoff = 0.05,#设定p值阈值
                           qvalueCutoff = 0.05,#设定q值阈值
                           readable = T)
    fourth_p_list[[paste0(antibody,"_",condition,"_fourth")]] <- barplot(fourth_GO,label_format = 50,showCategory = 10)+ggtitle(paste0(condition," fourth"))
    }
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_p <- plot_a_list(p_list,no_of_rows=2,no_of_cols=3)
# combined_p
ggsave("result/MEF_OE/histone_change_RNA_Cdkn2a.png",combined_p,width=15,height=10)


GO_p_list <- c(
  second_p_list[grepl("H3K27me3", names(second_p_list))],
  fourth_p_list[grepl("H3K27me3", names(fourth_p_list))]
)
combined_p <- plot_a_list(GO_p_list,no_of_rows=2,no_of_cols=3)
ggsave("result/MEF_OE/histone_change_RNA_GO_H3K27me3.png",combined_p,width=30,height=15)

GO_p_list <- c(
  second_p_list[grepl("H2AK119ub1", names(second_p_list))],
  fourth_p_list[grepl("H2AK119ub1", names(fourth_p_list))]
)
combined_p <- plot_a_list(GO_p_list,no_of_rows=2,no_of_cols=3)
ggsave("result/MEF_OE/histone_change_RNA_GO_H2AK119ub1.png",combined_p,width=30,height=15)

