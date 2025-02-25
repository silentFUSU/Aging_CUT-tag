rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(ChIPseeker)
library(EnsDb.Mmusculus.v79)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(ggalluvial)  
library(scales)
library(data.table)
# tissues <- c("mammarygland","BAT","CB","lung","kidney","aorta","brain","spleen",
             # "thymus","skin","bladder","bonemarrow","Hip","heart","muscle","iWAT")
tissues <- c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
             "thymus","skin","bladder","bonemarrow","Hip","heart",
             "muscle","jejunum","uterus","ovary","liver","tongue",
             "cecum","colon","testis","stomach","pancreas","iWAT","ileum")
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
summary <- data.frame()
for(tissue in tissues){
  region <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/10kb_all_significant_second_quadrant_after_remove_batch_effect.bed"))
  region <- region[which(region$V1 %in% paste0("chr",c(1:19))),]
  region$label <- paste0(region$V1,":",region$V2,"-",region$V3)
  region$tissue <- tissue_label_change(tissue)
  region <- region[,c("label","tissue")]
  summary <- rbind(summary,region)  
}

summary_count <- summary %>%   
  count(label)
summary_tissue <- summary %>%   
  group_by(label) %>%   
  summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
summary_count <- merge(summary_count,summary_tissue,by="label")
colnames(summary_count)[2] <- "count"
ggplot(summary_count, aes(x = count)) +  
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +  
  labs(x = "Count", y = "频数") +  
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.2) +
  ggtitle(paste0("H3K9me3 decrease and H3K27me3 increase common region")) +
  theme_minimal() +
  theme(  
    axis.title.x = element_text(size = 16),  # Increase x-axis label size  
    axis.title.y = element_text(size = 16),  # Increase y-axis label size  
    axis.text = element_text(size = 14),     # Increase axis tick labels size  
    plot.title = element_text(size = 18)     # Increase plot title size  
  )  


regions <- summary_count[c(summary_count$count >=7),]
bins <- read.table("~/ref_data/mm10_10kb_bins.bed")
bins$label <- paste0(bins$V1,":",bins$V2,"-",bins$V3)
bins <- bins[which(bins$label %in% regions$label),]
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
bins <- GRanges(seqnames = bins$V1,   
                     ranges = IRanges(start = bins$V2, end = bins$V3))
bins_anno <- annotatePeak(bins, tssRegion=c(-3000, 3000),
                               TxDb=txdb, annoDb="org.Mm.eg.db")
bins_anno <- unique(as.data.frame(bins_anno))
bins_anno <- bitr(bins_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
bins_anno_GO <- enrichGO( bins_anno$ENTREZID,
                              OrgDb = GO_database,
                              keyType = "ENTREZID",
                              ont = "BP",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05,
                              readable = T)

barplot(bins_anno_GO,label_format = 50)
to_plot <- data.frame()
for(tissue in tissues){
  H3K9me3_peak_region <- read.table(paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_10kb_in_young_merge-W1000-G3000-E100.bed"))
  H3K9me3_peak_region$label <- paste0(H3K9me3_peak_region$V1,":",H3K9me3_peak_region$V2,"-",H3K9me3_peak_region$V3)
  t_regions <- regions
  t_regions$condition <- "Out of H3K9me3 peaks"
  t_regions$condition[which(t_regions$label %in% H3K9me3_peak_region$label)] <- "In H3K9me3 peaks"
  t_to_plot <- as.data.frame(table(t_regions$condition))
  t_to_plot$percent <- t_to_plot$Freq/sum(t_to_plot$Freq) *100
  t_to_plot$tissue <- tissue_label_change(tissue)
  if(nrow(to_plot)==0){
    to_plot <- t_to_plot
  }else{
    to_plot <- rbind(to_plot,t_to_plot)
  }  
}
to_plot$position <- 100
to_plot$position[which(to_plot$Var1=="Out of H3K9me3 peaks")] <- to_plot[which(to_plot$Var1=="Out of H3K9me3 peaks"),"percent"]
to_plot_peaks <- to_plot[which(to_plot$Var1=="In H3K9me3 peaks"),]
to_plot_peaks <- to_plot_peaks[order(to_plot_peaks$percent),]

to_plot$tissue <- factor(to_plot$tissue,levels=to_plot_peaks$tissue)
total_label <- to_plot %>%  
  group_by(tissue) %>%  
  summarise(Freq_sum = sum(Freq))  
total_label$position <- 0
color <- setNames(c("#009980","#838B8B"),c("In H3K9me3 peaks","Out of H3K9me3 peaks"))
ggplot(to_plot, aes(x = tissue, y = percent, fill = Var1)) +  
  geom_bar(stat = 'identity',colour = "white") +   
  theme_minimal() +   
  scale_fill_manual(values = color) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20),legend.title = element_blank()) +
  geom_text(data = to_plot,   
            aes(label = Freq, y = position),   
            color = "black", size = 5, vjust = 0.5) +
  ylab("Proportion")+
  ggtitle("H3K27me3 increase and H3K9me3 decrease bins")

regions <- summary_count[c(summary_count$count >=7),]
bins <- read.table("~/ref_data/mm10_10kb_bins.bed")
bins$label <- paste0(bins$V1,":",bins$V2,"-",bins$V3)
bins <- bins[which(bins$label %in% regions$label),]
bins <- bins[,c("V1","V2","V3")]
bins$V2 <- bins$V2+1
setDT(bins)
setkey(bins,V1,V2,V3)
state_num <- "14"
p_list <- list()
for(tissue in tissues){
  file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_young[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_young <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_young <- chromHMM_young[which(chromHMM_young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_young$V2 <- chromHMM_young$V2+1
  chromHMM_young <- data.table(chromHMM_young)
  setDT(chromHMM_young) 
  setkey(chromHMM_young, V1, V2, V3) 
  
  file_dir <- paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/")  
  files_to_read <- list.files(path = file_dir, pattern = paste0(tissue, "_old[0-9]+_",state_num,"_segments_1k.bed"), full.names = TRUE)  
  file_list <- lapply(files_to_read,  read.delim, header = FALSE)  
  chromHMM_old <- Reduce(function(x, y) inner_join(x, y, by = c("V1", "V2", "V3", "V4")), file_list)  
  chromHMM_old <- chromHMM_old[which(chromHMM_old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  chromHMM_old$V2 <- chromHMM_old$V2+1
  chromHMM_old <- data.table(chromHMM_old)
  setDT(chromHMM_old) 
  setkey(chromHMM_old, V1, V2, V3) 
  
  overlaps_young <- foverlaps(bins, chromHMM_young, type = "any", nomatch = 0L)  
  overlaps_old <- foverlaps(bins, chromHMM_old, type = "any", nomatch = 0L)  
  
  overlaps_young$label <- paste0(overlaps_young$V1,":",overlaps_young$V2,"-",overlaps_young$V3)
  overlaps_old$label <- paste0(overlaps_old$V1,":",overlaps_old$V2,"-",overlaps_old$V3)
  overlaps <- merge(overlaps_young[,c("V4","label")],overlaps_old[,c("V4","label")],by="label")
  overlaps <- as.data.frame(overlaps)
  colnames(overlaps)[2:3] <- c("Young_state","Old_state") 
  to_plot <- overlaps %>%  
    group_by(Young_state, Old_state) %>%  
    summarize(freq = n())  
  to_plot$Young_state <- factor(to_plot$Young_state,levels=paste0("E",1:state_num))
  to_plot$Old_state <- factor(to_plot$Old_state,levels=paste0("E",1:state_num))
  colors <- hue_pal()(14) 
  colors <- setNames(colors,paste0("E",c(1:state_num)))
  
  p_list[[tissue]] <- ggplot(to_plot, aes(axis1 = Young_state, axis2 = Old_state, y = freq)) +  
    geom_alluvium(aes(fill = Young_state)) +  
    scale_fill_manual(values = colors) +
    geom_stratum() +  
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
    theme_minimal() +  
    labs(y = "Count", x = "State Transition", 
         fill = "Young State")+
    ggtitle(paste0(tissue_label_change(tissue)))
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_plot <- plot_a_list(p_list,no_of_cols = 7,no_of_rows = ceiling(length(p_list)/7))
ggsave(paste0("result/all/H3K27me3_H3K9me3/all_tissues_second_quadrant_after_remove_batch_effect_common_region_without_chrXY_larger10_chromHMM_annotation.png"),width = 35,height = ceiling(length(p_list)/7)*5, type="cairo")


regions <- summary_count[c(summary_count$count >=10),]
bins <- read.table("~/ref_data/mm10_10kb_bins.bed")
bins$label <- paste0(bins$V1,":",bins$V2,"-",bins$V3)
bins <- bins[which(bins$label %in% regions$label),]
H3K9me3_peaks <- read.table("data/samples/BAT/H3K9me3/bed/H3K9me3_young_merge-W1000-G3000-E100_compress.bed")

bins$V2 <- bins$V2+1
H3K9me3_peaks$V2 <- H3K9me3_peaks$V2 +1

bins <- as.data.table(bins)
H3K9me3_peaks <- as.data.table(H3K9me3_peaks)
setDT(bins)
setkey(bins,V1,V2,V3)

setDT(H3K9me3_peaks)
setkey(H3K9me3_peaks,V1,V2,V3)

overlaps <- foverlaps(bins, H3K9me3_peaks, type = "any", nomatch = 0L)  
overlaps <- overlaps[,c("V1","V2","V3","V4")]
overlaps <- overlaps[!duplicated(overlaps$V4), ] 
