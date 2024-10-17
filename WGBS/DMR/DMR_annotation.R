rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(genomation)
library(methylKit)
library(ChIPseeker)
library(ggplot2)
library(stringr)
library(dplyr)

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}

cpg.file="~/ref_data/for_normal_mapping/mm10/cpgi.mm10.bed.txt"
# cpg <- read.table(cpg.file)
# write.table(cpg,"~/ref_data/for_normal_mapping/mm10/cpgi.mm10.bed.txt",sep="\t",append = F,quote = F,row.names = F,col.names = F)
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
    }
  }
  return(tissue_label)
}

DMR_annotation <- function(tissue){
  p_list <- list()
  DMR <- read.table(paste0("data/samples/WGBS/",tissue,"/DSS_table/",tissue,"_DMR_delta0.txt"),header = T)
  DMR <- DMR[which(DMR$chr %in% paste0("chr",c(1:19,"X","Y"))),]
  cpg.shore.obj=readFeatureFlank(cpg.file,flank = 2000,feature.flank.name=c("CpGi","shores"))
  cpg.shelf.obj=readFeatureFlank(cpg.file,flank = 4000,feature.flank.name=c("CpGi","shelf"))
  colnames(DMR)[1:3] <- c("chr","start","end")
  increase <- DMR[which(DMR$diff.Methy > 0),c(1:3)] 
  decrease <- DMR[which(DMR$diff.Methy < 0),c(1:3)]
  
  shore_increase_anno=annotateWithFeatureFlank(as(increase,"GRanges"),
                                              cpg.shore.obj$CpGi,cpg.shore.obj$shores,
                                              feature.name="CpGi",flank.name="shores")
  shore_decrease_anno=annotateWithFeatureFlank(as(decrease,"GRanges"),
                                               cpg.shore.obj$CpGi,cpg.shore.obj$shores,
                                               feature.name="CpGi",flank.name="shores")
  shelf_increase_anno=annotateWithFeatureFlank(as(increase,"GRanges"),
                                               cpg.shelf.obj$CpGi,cpg.shelf.obj$shelf,
                                               feature.name="CpGi",flank.name="shelf")
  shelf_decrease_anno=annotateWithFeatureFlank(as(decrease,"GRanges"),
                                               cpg.shelf.obj$CpGi,cpg.shelf.obj$shelf,
                                               feature.name="CpGi",flank.name="shelf")
  increase_summary <- data.frame(location = c("Island","Shore","Shelf","Open sea"),
                                 percent = c(shore_increase_anno@precedence[[1]],
                                             shore_increase_anno@precedence[[2]],
                                             shelf_increase_anno@precedence[[2]]-shore_increase_anno@precedence[[2]],
                                             shelf_increase_anno@precedence[[3]]))
  increase_summary$condition <- "Hyper"
  
  decrease_summary <- data.frame(location = c("Island","Shore","Shelf","Open sea"),
                                 percent = c(shore_decrease_anno@precedence[[1]],
                                             shore_decrease_anno@precedence[[2]],
                                             shelf_decrease_anno@precedence[[2]]-shore_decrease_anno@precedence[[2]],
                                             shelf_decrease_anno@precedence[[3]]))
  decrease_summary$condition <- "Hypo"
  
  summary <- rbind(increase_summary,decrease_summary)
  summary$location <- factor(summary$location, levels = c("Open sea","Shelf","Shore","Island"))
  color <- setNames(c("#009980","#739940","#E69900","#838B8B"),c("Island","Shore","Shelf","Open sea"))
  p_list[[1]] <- ggplot( summary, aes(x = condition, y = percent, fill = location)) +  
    geom_bar(stat = 'identity',color="white") +   
    theme_minimal() +   
    scale_fill_manual(values = color) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue)))
  
  GO_database <- 'org.Mm.eg.db'
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  increase_Granges <- GRanges(seqnames =increase$chr,   
          ranges = IRanges(start = increase$start, end = increase$end))
  increase_anno <- annotatePeak(increase_Granges, tssRegion=c(-3000, 3000),
                               TxDb=txdb, annoDb="org.Mm.eg.db")
  
  decrease_Granges <- GRanges(seqnames =decrease$chr,   
                              ranges = IRanges(start = decrease$start, end = decrease$end))
  decrease_anno <- annotatePeak(decrease_Granges, tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Mm.eg.db") 
  increase_anno_summary <- increase_anno@annoStat
  increase_anno_summary$condition <- "Hyper"
  exon_sum <- increase_anno_summary %>%
    filter(grepl("Exon",Feature)) %>%
    summarise(Frequency = sum(Frequency)) 
  new_exon_row <- data.frame(Feature = "Exon", Frequency = exon_sum$Frequency, condition="Hyper")  
  
  intron_sum <- increase_anno_summary %>%
    filter(grepl("Intron",Feature)) %>%
    summarise(Frequency = sum(Frequency)) 
  new_intron_row <- data.frame(Feature = "Intron", Frequency = intron_sum$Frequency, condition="Hyper")  
  
  prom_sum <- increase_anno_summary %>%
    filter(grepl("Promoter",Feature)) %>%
    summarise(Frequency = sum(Frequency)) 
  Distal_prom <-  data.frame(Feature = "Distal prom", Frequency = (prom_sum$Frequency - increase_anno_summary$Frequency[which(increase_anno_summary$Feature == "Promoter (<=1kb)")]), condition="Hyper")  
  
  increase_anno_summary <- rbind(increase_anno_summary,new_exon_row)
  increase_anno_summary <- rbind(increase_anno_summary,new_intron_row)
  increase_anno_summary <- rbind(increase_anno_summary,Distal_prom)
  
  decrease_anno_summary <- decrease_anno@annoStat
  decrease_anno_summary$condition <- "Hypo"
  exon_sum <- decrease_anno_summary %>%
    filter(grepl("Exon",Feature)) %>%
    summarise(Frequency = sum(Frequency)) 
  new_exon_row <- data.frame(Feature = "Exon", Frequency = exon_sum$Frequency, condition="Hyper")  
  
  intron_sum <- decrease_anno_summary %>%
    filter(grepl("Intron",Feature)) %>%
    summarise(Frequency = sum(Frequency)) 
  new_intron_row <- data.frame(Feature = "Intron", Frequency = intron_sum$Frequency, condition="Hyper")  
  
  prom_sum <- decrease_anno_summary %>%
    filter(grepl("Promoter",Feature)) %>%
    summarise(Frequency = sum(Frequency)) 
  Distal_prom <-  data.frame(Feature = "Distal prom", Frequency = (prom_sum$Frequency - decrease_anno_summary$Frequency[which(decrease_anno_summary$Feature == "Promoter (<=1kb)")]), condition="Hyper")  
  
  decrease_anno_summary <- rbind(decrease_anno_summary,new_exon_row)
  decrease_anno_summary <- rbind(decrease_anno_summary,new_intron_row)
  decrease_anno_summary <- rbind(decrease_anno_summary,Distal_prom)
  
  colnames(increase_anno_summary)[2] <- "Hyper"
  colnames(decrease_anno_summary)[2] <- "Hypo"
  anno_summary <- merge(increase_anno_summary[which(increase_anno_summary$Feature %in% c("Distal Intergenic","Downstream (<=300)","3' UTR","Intron","Exon","5' UTR","Promoter (<=1kb)","Distal prom")),c(1:2)],
                        decrease_anno_summary[which(decrease_anno_summary$Feature %in% c("Distal Intergenic","Downstream (<=300)","3' UTR","Intron","Exon","5' UTR","Promoter (<=1kb)","Distal prom")),c(1:2)],
                        by="Feature")
  
  anno_summary$Feature <- as.character(anno_summary$Feature)
  anno_summary$Feature[which(anno_summary$Feature=="Downstream (<=300)")] <- "Downstream"
  anno_summary$Feature[which(anno_summary$Feature=="Distal Intergenic")] <- "Intergenic"
  anno_summary$Feature[which(anno_summary$Feature=="Promoter (<=1kb)")] <- "Promoter"
  anno_summary$Feature <- factor(anno_summary$Feature, levels = c("Intergenic", "Downstream", "3' UTR", "Intron", "Exon", "5' UTR", "Promoter", "Distal prom"))
  color <- setNames(c("#838B8B","#E69900","#BF9915","#99992A","#739940","#4C9955","#26996A","#009980"),c("Intergenic", "Downstream", "3' UTR", "Intron", "Exon", "5' UTR", "Promoter", "Distal prom"))
  anno_summary_to_plot <- reshape2::melt(anno_summary)
  p_list[[2]] <- ggplot( anno_summary_to_plot, aes(x = variable, y = value, fill = Feature)) +  
    geom_bar(stat = 'identity',colour = "white") +   
    theme_minimal() +   
    scale_fill_manual(values = color) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue)))
  return(p_list)
  }

CpG_annotation <- list()
gene_annotation <- list()

tissues <-sort(c("liver","lung","kidney","ileum","Hip","mammarygland","skin","bonemarrow","jejunum","colon","ovary","CB","BAT","thymus","testis"))
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  p_list <- DMR_annotation(tissues[i])
  CpG_annotation[[i]] <- p_list[[1]]
  names(CpG_annotation)[i] <- tissue_label_change(tissue)
  gene_annotation[[i]] <- p_list[[2]]
  names(gene_annotation)[i] <- tissue_label_change(tissue)
}
CpG_annotation_sort <- CpG_annotation[sort$tissue]

CpG_annotation_combine <- plot_a_list(CpG_annotation_sort,no_of_rows = 3,no_of_cols = 5) + patchwork::plot_annotation(title = "CpG context",theme = theme(plot.title = element_text(size = 40,hjust = 0.5)))  
ggsave("result/WGBS/all_tissues_DMR_delta0_annotation_CpG_context.png",CpG_annotation_combine,width = 15,height = 15,type="cairo")

gene_annotation_sort <- gene_annotation[sort$tissue]
gene_annotation_combine <- plot_a_list(gene_annotation_sort,no_of_rows = 3,no_of_cols = 5)+ patchwork::plot_annotation(title = "Gene location",theme = theme(plot.title = element_text(size = 40,hjust = 0.5)))  
ggsave("result/WGBS/all_tissues_DMR_delta0_annotation_gene_location.png",gene_annotation_combine,width = 15,height = 15,type="cairo")






