rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(reshape2)
library(ChIPseeker)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(stringr)
library(data.table)
library(rtracklayer)  
library(ggrepel)  
library(GenomicFeatures)  
tissue <- "lung"
antibody <- "H3K27me3"
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
bin_size <- function(antibody){
  if(antibody %in% c("H3K9me3","H3K27me3","H3K36me3")){
    return("10kb")
  }else{
    return("1kb")
  }
}
rna_name <- function(tissue){
  return(tissue)
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
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
GO_database <- 'org.Mm.eg.db'

txdb <- makeTxDbFromGFF("/storage/zhangyanxiaoLab/share/gtf/mm10.gencode.vM25.annotation.gtf",format = "gtf")
gtf_data <- import("/storage/zhangyanxiaoLab/share/gtf/mm10.gencode.vM25.annotation.gtf", format = "gtf")
gene_table <- as.data.frame(gtf_data@elementMetadata[,c("transcript_id","gene_name")])
gene_table <- gene_table[!duplicated(gene_table$transcript_id), ]
colnames(gene_table) <- c("transcriptId","SYMBOL")
Histone_relationship_with_RNA <- function(tissue,antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  histone <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  histone <- histone[,c("Chr","Start","End","LogFC.old.young","Significant")]
  peak_obj <- GRanges(seqnames = histone$Chr,   
                      ranges = IRanges(start = histone$Start, end = histone$End))
  peak_anno <- annotatePeak(peak_obj, tssRegion=c(-1000, 1000),
                            TxDb=txdb)
  peak_anno <- unique(as.data.frame(peak_anno))
  peak_anno <- peak_anno[which(str_detect(peak_anno$annotation,"Promoter")),]
  peak_anno <- peak_anno[which(peak_anno$distanceToTSS==0),]
  peak_anno <- merge(peak_anno,gene_table,by="transcriptId")
  peak_anno <- peak_anno %>%  
    group_by(SYMBOL) %>%  
    slice_max(order_by = geneLength, n = 1, with_ties = FALSE) %>%
    ungroup()  
  
  rna <- read.csv(paste0("data/samples/RNA/",rna_name(tissue),"/diff_expression_gene.csv"))
  peak_anno$label <- paste0(peak_anno$seqnames,"-",peak_anno$start,"-",peak_anno$end)
  rna <- merge(rna[,c("X","logFC","Significant")],peak_anno[,c("SYMBOL","label")],by.x="X",by.y="SYMBOL")
  colnames(rna) <- c("Gene","RNA_logFC","RNA_Significant","label")
  histone$label <- paste0(histone$Chr,"-",histone$Start,"-",histone$End)
  to_plot <- merge(rna,histone[,c("label","LogFC.old.young","Significant")],by="label")
  # to_plot <- to_plot[which(to_plot$Significant!="Stable"),]
  if(nrow(to_plot)==0){
    return(0)
  }
  to_plot$tissue <- tissue_label_change(tissue)
  # gene_summary <<- rbind(gene_summary,to_plot[which(to_plot$LogFC <0 & to_plot$RNA_Significant=="Up"),c("Gene","tissue")])
  top_genes <- to_plot[which(to_plot$RNA_Significant=="Up"),] %>%  
    arrange(desc(RNA_logFC)) %>%  
    head(5)  
  bottom_genes <- to_plot[which(to_plot$RNA_Significant=="Down"),] %>%  
    arrange(RNA_logFC) %>%  
    head(5)  
  highlight_genes <- rbind(top_genes, bottom_genes)  
  x_range <- range(to_plot$RNA_logFC, na.rm = TRUE)  
  y_range <- range(to_plot$LogFC.old.young, na.rm = TRUE)  
  x_pos_right <- x_range[2] * 0.9    
  x_pos_left <- x_range[1] * 0.9   
  y_pos_top <- y_range[2] * 0.9    
  y_pos_bottom <- y_range[1] * 0.9 
  ggplot() +
    geom_point(data=to_plot, mapping=aes(RNA_logFC, LogFC.old.young),color = "grey",alpha=0.5) +  
    geom_point(data=to_plot[which(to_plot$RNA_Significant=="Up"),], mapping=aes(RNA_logFC, LogFC.old.young),color = "red") +
    geom_point(data=to_plot[which(to_plot$RNA_Significant=="Down"),], mapping=aes(RNA_logFC, LogFC.old.young),color = "blue") +
    # geom_point(data=to_plot[which(to_plot$RNA_logFC>0 & to_plot$LogFC<0 & to_plot$RNA_Significant=="Up"),], mapping=aes(RNA_logFC, LogFC),color = "#48466d") +
    # geom_point(data=to_plot[which(to_plot$RNA_logFC>0 & to_plot$LogFC>0),], mapping=aes(RNA_logFC, LogFC),color = "#00b8a9") +
    # geom_point(data=to_plot[which(to_plot$RNA_logFC<0 & to_plot$LogFC<0),], mapping=aes(RNA_logFC, LogFC),color = "#ff9a00") +
    # 坐标轴
    labs(x="RNA log2(Fold Change)",
         y=paste0(antibody," log2(Fold Change)")) +
    geom_point(color="grey")+
    theme(legend.position = "bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme_bw()+theme(text = element_text(size = 18))+
    ggtitle(tissue_label_change(tissue),paste0(antibody," all regions"))+
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_text_repel(data = highlight_genes, aes(RNA_logFC, LogFC.old.young, label = Gene),   
                    size = 5, # 字体大小  
                    nudge_y = 0.2)+
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$RNA_logFC > 0 & to_plot$LogFC.old.young > 0 & to_plot$RNA_Significant == "Up"), ])),  
             x = x_pos_right, y = y_pos_top, colour = "#00b8a9", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$RNA_logFC < 0 & to_plot$LogFC.old.young < 0 & to_plot$RNA_Significant == "Down"), ])),  
             x = x_pos_left, y = y_pos_bottom, colour = "#ff9a00", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$RNA_logFC < 0 & to_plot$LogFC.old.young > 0 & to_plot$RNA_Significant == "Down"), ])),  
             x = x_pos_left, y = y_pos_top, colour = "#f6416c", size = 5) +  
    annotate("text", label = paste0(nrow(to_plot[which(to_plot$RNA_logFC > 0 & to_plot$LogFC.old.young < 0 & to_plot$RNA_Significant == "Up"), ])),  
             x = x_pos_right, y = y_pos_bottom, colour = "#48466d", size = 5)  
  }
tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K4me3","H3K4me1","H3K27ac")
# antibodys <- c("H3K4me3","H3K4me1","H3K27ac")
for(antibody in antibodys){
  p_list <- list()
  # gene_summary_increase <- data.frame()
  for(tissue in tissues){
    p_list[[tissue]] <- Histone_relationship_with_RNA(tissue,antibody)
  }
  p_list <- p_list[sapply(p_list, function(x) !is.null(x) && !is.numeric(x) || (is.numeric(x) && x != 0))]   
  # gene_summary_count <- gene_summary %>%
  #   count(Gene)
  # gene_summary_tissue <- gene_summary %>%   
  #   group_by(Gene) %>%   
  #   summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
  # 
  # gene_summary_count <- merge(gene_summary_tissue,gene_summary_count,by="Gene")
  # 
  combined_plot <- plot_a_list(p_list,no_of_rows = 4,no_of_cols = 7)
  ggsave(paste0("result/RNA/histone_relationship_with_RNA/",antibody,"/",antibody,"_relationship_with_RNA_bin_level_all_regions.png"),combined_plot,width = 34,height = 20,type="cairo")
}


# txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene

Histone_relationship_with_RNA_peak_level <- function(tissue,antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  histone <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff_after_remove_batch_effect.csv"))
  peaks <- read.table(paste0("data/samples/",tissue,"/",antibody,"/bed/",antibody,"_young_old_merge-W1000-G3000-E100_compress.bed"))
  peaks$V2 <- peaks$V2+1
  histone$Start <- histone$Start + 1
  histone <- histone[,c("Chr","Start","End","LogFC.old.young","Significant")]
  histone <- as.data.table(histone)
  peaks <- as.data.table(peaks)
  setDT(histone)
  setDT(peaks)
  setkey(histone,Chr,Start,End)
  setkey(peaks,V1,V2,V3)
  overlap <- foverlaps(histone, peaks, type = "any", nomatch = 0L)  
  result <- overlap %>%  
    group_by(V4) %>%  
    summarise(  
      total = n(),  
      up_count = sum(Significant == "Up"),  
      down_count = sum(Significant == "Down"),  
      up_ratio = up_count / total,  
      down_ratio = down_count / total,
      LogFC = median(LogFC.old.young, na.rm = TRUE)
    ) %>%  
    mutate(  
      Classification = case_when(  
        up_ratio > 0.3 & LogFC>0 ~ "increase",  
        down_ratio > 0.3 & LogFC<0~ "decrease",  
        TRUE ~ "stable"  
      )  
    ) %>%  
    select(V4, Classification,LogFC)  
  peaks <- merge(peaks,result,by="V4")
  peak_obj <- GRanges(seqnames = peaks$V1,   
                      ranges = IRanges(start = peaks$V2, end = peaks$V3))
  peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                            TxDb=txdb)
  peak_anno <- unique(as.data.frame(peak_anno))
  peak_anno <- peak_anno[which(str_detect(peak_anno$annotation,"Promoter")),]
  peak_anno <- peak_anno[which(peak_anno$distanceToTSS==0),]
  peak_anno <- merge(peak_anno,gene_table,by="transcriptId")
  peak_anno <- peak_anno %>%  
    group_by(SYMBOL) %>%  
    slice_max(order_by = geneLength, n = 1, with_ties = FALSE) %>%
    ungroup()  
  
  rna <- read.csv(paste0("data/samples/RNA/",rna_name(tissue),"/diff_expression_gene.csv"))
  peak_anno$label <- paste0(peak_anno$seqnames,"-",peak_anno$start,"-",peak_anno$end)
  rna <- merge(rna[,c("X","logFC","Significant")],peak_anno[,c("SYMBOL","label")],by.x="X",by.y="SYMBOL")
  colnames(rna) <- c("Gene","RNA_logFC","RNA_Significant","label")
  peaks$label <- paste0(peaks$V1,"-",peaks$V2,"-",peaks$V3)
  to_plot <- merge(rna,peaks[,c("label","LogFC","Classification")],by="label")
  to_plot <- to_plot[which(to_plot$Classification!="stable"),]
  to_plot$tissue <- tissue_label_change(tissue)
  gene_summary <<- rbind(gene_summary,to_plot[which(to_plot$LogFC <0 & to_plot$RNA_Significant=="Up"),c("Gene","tissue")])
  top_genes <- to_plot[which(to_plot$RNA_Significant=="Up"),] %>%  
    arrange(desc(RNA_logFC)) %>%  
    head(5)  
  bottom_genes <- to_plot[which(to_plot$RNA_Significant=="Down"),] %>%  
    arrange(RNA_logFC) %>%  
    head(5)  
  highlight_genes <- rbind(top_genes, bottom_genes)  
  
  p <- ggplot() +
    geom_point(data=to_plot, mapping=aes(RNA_logFC, LogFC),color = "grey",alpha=0.5) +  
    geom_point(data=to_plot[which(to_plot$RNA_Significant=="Up"),], mapping=aes(RNA_logFC, LogFC),color = "red") +
    geom_point(data=to_plot[which(to_plot$RNA_Significant=="Down"),], mapping=aes(RNA_logFC, LogFC),color = "blue") +
    # geom_point(data=to_plot[which(to_plot$RNA_logFC>0 & to_plot$LogFC<0 & to_plot$RNA_Significant=="Up"),], mapping=aes(RNA_logFC, LogFC),color = "#48466d") +
    # geom_point(data=to_plot[which(to_plot$RNA_logFC>0 & to_plot$LogFC>0),], mapping=aes(RNA_logFC, LogFC),color = "#00b8a9") +
    # geom_point(data=to_plot[which(to_plot$RNA_logFC<0 & to_plot$LogFC<0),], mapping=aes(RNA_logFC, LogFC),color = "#ff9a00") +
    # 坐标轴
    labs(x="RNA log2(Fold Change)",
         y=paste0(antibody," log2(Fold Change)")) +
    geom_point(color="grey")+
    theme(legend.position = "bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme_bw()+theme(text = element_text(size = 18))+
    ggtitle(tissue_label_change(tissue),paste0(antibody," unstable regions"))+
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_text_repel(data = highlight_genes, aes(RNA_logFC, LogFC, label = Gene),   
                  size = 5, # 字体大小  
                  nudge_y = 0.2)
    # annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young>0 & peak_anno$logFC>0),])),x=10, y=4,colour="#00b8a9",size=5)+
    # annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young<0 & peak_anno$logFC<0),])),x=-5, y=-4,colour="#ff9a00",size=5)+
    # annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young>0 & peak_anno$logFC<0),])),x=-5, y=4,colour="#f6416c",size=5)+
    # annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young<0 & peak_anno$logFC>0),])),x=10, y=-4,colour="#48466d",size=5)
  return(p)
  }


tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
p_list <- list()
antibody <- "H3K27me3"
gene_summary <- data.frame()
for(tissue in tissues){
  p_list[[tissue]] <- Histone_relationship_with_RNA_peak_level(tissue,antibody)
}

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_plot <- plot_a_list(p_list,no_of_rows = 4,no_of_cols = 7)
ggsave(paste0("result/RNA/histone_relationship_with_RNA/",antibody,"/",antibody,"_relationship_with_RNA_peak_level.png"),combined_plot,width = 42,height = 24,type="cairo")

gene_summary_count <- gene_summary %>%   
  count(Gene)
gene_summary_tissue <- gene_summary %>%   
  group_by(Gene) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  

gene_summary_count <- merge(gene_summary_tissue,gene_summary_count,by="Gene")
