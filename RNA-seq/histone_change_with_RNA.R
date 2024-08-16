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
tissue <- "CB"
antibody <- "H3K27me3"
bin_size <- function(antibody){
  if(antibody %in% c("H3K9me3","H3K27me3","H3K36me3")){
    return("10kb")
  }else{
    return("1kb")
  }
}
rna_name <- function(tissue){
  if(tissue=="brain"){
    return("FC")
  }
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
    }
  }
  return(tissue_label)
} 
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}

Histone_relationship_with_RNA <- function(antibody,tissue){
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  GO_database <- 'org.Mm.eg.db'
  if(antibody=="ATAC"){
    data_path <- "data/samples/ATAC/"
  }else{
    data_path <- "data/samples/"
  }
  histone <- read.csv(paste0(data_path,tissue,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff_after_remove_batch_effect.csv"))
  rna <- read.csv(paste0("data/samples/RNA/",rna_name(tissue),"/diff_expression_gene_nodup.csv"))
  
  rna_gene_change <- rna$X[which(rna$Significant != "Stable")]
  colnames(rna)[1]<-"SYMBOL"
  peak_obj <- GRanges(seqnames = histone$Chr,   
                      ranges = IRanges(start = histone$Start, end = histone$End))
  peak_anno <- annotatePeak(peak_obj, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  peak_anno <- unique(as.data.frame(peak_anno))
  peak_anno <- peak_anno[which(str_detect(peak_anno$annotation,"Promoter")),]
  peak_anno <- peak_anno[which(peak_anno$SYMBOL %in% rna_gene_change),]
  
  peak_anno$label <- paste0(peak_anno$seqnames,"-",peak_anno$start,"-",peak_anno$end)
  histone$label <- paste0(histone$Chr,"-",histone$Start,"-",histone$End)
  peak_anno <- merge(peak_anno,histone[,c("LogFC.old.young","label","Significant_bar")],by="label")
  peak_anno <- merge(peak_anno,rna[,c("SYMBOL","logFC")],by="SYMBOL")
  peak_anno <- peak_anno[which(peak_anno$Significant_bar != "Stable" & peak_anno$distanceToTSS==0),]
  t_sort_table <- data.frame(tissue = tissue,count = nrow(peak_anno[which(peak_anno$LogFC.old.young>0 & peak_anno$logFC>0),]))
  sort_table <<- rbind(sort_table,t_sort_table)
  p<- ggplot() +
    geom_point(data=peak_anno, mapping=aes(logFC, LogFC.old.young),color = "grey",alpha=0.5) +  
    # geom_point(data=peak_anno[which(peak_anno$LogFC.old.young<0 & peak_anno$logFC>0),], mapping=aes(logFC, LogFC.old.young),color = "#48466d") +
    geom_point(data=peak_anno[which(peak_anno$LogFC.old.young>0 & peak_anno$logFC>0),], mapping=aes(logFC, LogFC.old.young),color = "#00b8a9") +
    geom_point(data=peak_anno[which(peak_anno$LogFC.old.young<0 & peak_anno$logFC<0),], mapping=aes(logFC, LogFC.old.young),color = "#ff9a00") +
    # 坐标轴
    labs(x="RNA log2(Fold Change)",
         y=paste0(antibody," log2(Fold Change)")) +
    geom_point(color="grey")+
    theme(legend.position = "bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme_bw()+theme(text = element_text(size = 18))+
    ggtitle(tissue_label_change(tissue))+
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young>0 & peak_anno$logFC>0),])),x=10, y=4,colour="#00b8a9",size=5)+
    annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young<0 & peak_anno$logFC<0),])),x=-5, y=-4,colour="#ff9a00",size=5)+
    annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young>0 & peak_anno$logFC<0),])),x=-5, y=4,colour="#f6416c",size=5)+
    annotate("text",label = paste0(nrow(peak_anno[which(peak_anno$LogFC.old.young<0 & peak_anno$logFC>0),])),x=10, y=-4,colour="#48466d",size=5)
  print(p)
  return(p)
}

tissues <-c("thymus","CB","uterus","lung","muscle","skin","spleen","bonemarrow","heart","liver","kidney","testis","Hip","brain", "ileum","aorta","tongue","bladder","stomach","jejunum","colon","cecum")

antibody <- "ATAC"
p_list <- list()
sort_table <- data.frame(tissue = as.character(),
                         count = as.numeric())

for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  p_list[[i]] <- Histone_relationship_with_RNA(antibody,tissue)
  names(p_list)[i] <- tissue
}
sort_table_order <- sort_table[order(-sort_table$count),]
p_list_order <- p_list[sort_table_order$tissue]
combined_plot <- plot_a_list(p_list_order, 4, 6)
ggsave(paste0("result/RNA/histone_relationship_with_RNA/",antibody,"/",antibody,"_relationship_with_RNA.png"),combined_plot,width = 30,height = 20,type="cairo")

p_list <- list()
antibodys <- c("H3K27ac","H3K4me3","ATAC")
for(i in c(1:length(antibodys))){
  antibody <- antibodys[i]
  p_list[[i]] <- Histone_relationship_with_RNA(antibody,"CB")
}
combined_plot <- plot_a_list(p_list, 1,3)
print(combined_plot)
