rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(ggalluvial)
library(GenomicRanges)
tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver","ileum",
                  "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT"))
state_num <- 15
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
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

summary <- data.frame()
for(tissue in tissues){
  young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
  young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
  young1$label <- paste(young1$V1,young1$V2,young1$V3,young1$V4,sep = "-")
  young2$label <- paste(young2$V1,young2$V2,young2$V3,young2$V4,sep = "-")
  young <- young1[which(young1$label %in% intersect(young1$label,young2$label)),]
  
  old1 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_old1_",state_num,"_segments_1k.bed"),header = F)
  old2 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_old2_",state_num,"_segments_1k.bed"),header = F)
  old1$label <- paste(old1$V1,old1$V2,old1$V3,old1$V4,sep = "-")
  old2$label <- paste(old2$V1,old2$V2,old2$V3,old2$V4,sep = "-")
  old <- old1[which(old1$label %in% intersect(old1$label,old2$label)),]
  
  young <- as.data.frame(table(young$V4))
  old <- as.data.frame(table(old$V4))
  t_summary <- merge(young,old,by="Var1")
  t_summary$mean <- rowMeans(t_summary[,-1])
  t_summary$percentage <- t_summary$mean/sum(t_summary$mean) *100
  t_summary$coverage <- t_summary$mean * 1000
  t_summary$tissue <- tissue_label_change(tissue)
  t_summary <- t_summary[,c("Var1","percentage","coverage","tissue")]
  summary <- rbind(summary,t_summary)
}

to_plot <- summary
dictionary <- list("E1"=1, "E2"=2, "E3"=3,
                   "E4"=4, "E5"=5, "E6"=7,
                   "E7"=8, "E8"=6, "E9"=9,
                   "E10"=10,"E11"=11,"E12"=15,
                   "E13"=12,"E14"=13,"E15"=14)
keys <- names(dictionary)
values <- unlist(dictionary)
to_plot$Var1 <- values[match(to_plot$Var1, keys)]
to_plot$Var1 <- paste0("E",to_plot$Var1)
to_plot$Var1 <- gsub("^E", "state", to_plot$Var1)
to_plot$Var1 <- factor(to_plot$Var1,rev(paste0("state",1:state_num)))
color <- read.table("data/samples/20_distinct_color.txt")
color <- setNames(color$V1,paste0("state",1:state_num))
ggplot(to_plot, aes(x = percentage, y = Var1)) +
  geom_boxplot(fill="gray") +
  theme_bw()+
  # scale_fill_manual(values=color,breaks = paste0("state",1:state_num))+
  labs(x = "Percentage", y = NULL) +
  theme_minimal()
breaks = seq(6,9,1)
labels <- ifelse(10**breaks >= 1e9, paste0(round(10**breaks/1e9, 1), "Gb"),
                 ifelse(10**breaks >= 1e6, paste0(round(10**breaks/1e6, 1), "Mb"), 
                        ifelse(10**breaks >= 1e3, paste0(round(10**breaks/1e3, 1), "Kb"),
                               as.character(10**breaks))))
p <- ggplot(to_plot, aes(x = log10(coverage), y = Var1),) +
  geom_boxplot(fill="gray") +
  # scale_fill_manual(values=color,breaks = paste0("state",1:state_num))+
  scale_x_continuous(breaks=breaks,
                     labels=labels )+
  labs(x = "Coverage", y = NULL) +
  theme_bw()

ggsave("result/figures/chromHMM_15_state_coverage.pdf",p,width = 6,height = 10)
### annotation
library(org.Mm.eg.db)
library(ChIPseeker)
library("AnnotationDbi")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
young_summary <- data.frame()
old_summary <- data.frame()
for(tissue in tissues){
  young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
  young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
  young1$label <- paste(young1$V1,young1$V2,young1$V3,young1$V4,sep = "-")
  young2$label <- paste(young2$V1,young2$V2,young2$V3,young2$V4,sep = "-")
  young <- young1[which(young1$label %in% intersect(young1$label,young2$label)),]
  young <- young[which(young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  young$label <- paste(young$V1,young$V2,young$V3,sep = "-")
  young_bin <- GRanges(seqnames = young$V1,   
                  ranges = IRanges(start = young$V2, end = young$V3))  
  
  young_bin_anno <- annotatePeak(young_bin, tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Mm.eg.db")
  young_bin_anno <- as.data.frame(young_bin_anno@anno)
  young_bin_anno$annotation_label <- ifelse(grepl("Exon|Intron", young_bin_anno$annotation),
                                      gsub("^(Exon|Intron).*", "\\1", young_bin_anno$annotation),
                                      young_bin_anno$annotation)      
  young_bin_anno$label <- paste(young_bin_anno$seqnames,young_bin_anno$start,young_bin_anno$end,sep="-")
  young_bin_anno <- young_bin_anno[,c("label","annotation_label")]
  young_bin_anno <- merge(young_bin_anno, young[,c("label","V4")], by="label")
  young_result <- young_bin_anno %>%
    group_by(V4, annotation_label) %>%
    summarise(count = n()) %>%
    mutate(percentage = (count / sum(count)) * 100)
  young_result$tissue <- tissue_label_change(tissue)
  young_summary <- rbind(young_summary,young_result)
  
  old1 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_old1_",state_num,"_segments_1k.bed"),header = F)
  old2 <- read.delim(paste0("result/all/ChromHMM/all_tissues_normal_chr/",state_num,"_all_tissues/split_1k/",tissue,"_old2_",state_num,"_segments_1k.bed"),header = F)
  old1$label <- paste(old1$V1,old1$V2,old1$V3,old1$V4,sep = "-")
  old2$label <- paste(old2$V1,old2$V2,old2$V3,old2$V4,sep = "-")
  old <- old1[which(old1$label %in% intersect(old1$label,old2$label)),]
  old <- old[which(old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
  old$label <- paste(old$V1,old$V2,old$V3,sep = "-")
  old_bin <- GRanges(seqnames = old$V1,   
                       ranges = IRanges(start = old$V2, end = old$V3))  
  
  old_bin_anno <- annotatePeak(old_bin, tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Mm.eg.db")
  old_bin_anno <- as.data.frame(old_bin_anno@anno)
  old_bin_anno$annotation_label <- ifelse(grepl("Exon|Intron", old_bin_anno$annotation),
                                            gsub("^(Exon|Intron).*", "\\1", old_bin_anno$annotation),
                                            old_bin_anno$annotation)      
  old_bin_anno$label <- paste(old_bin_anno$seqnames,old_bin_anno$start,old_bin_anno$end,sep="-")
  old_bin_anno <- old_bin_anno[,c("label","annotation_label")]
  old_bin_anno <- merge(old_bin_anno, old[,c("label","V4")], by="label")
  old_result <- old_bin_anno %>%
    group_by(V4, annotation_label) %>%
    summarise(count = n()) %>%
    mutate(percentage = (count / sum(count)) * 100)
  old_result$tissue <- tissue_label_change(tissue)
  old_summary <- rbind(old_summary,old_result)
  }

to_plot <- as.data.frame(old_summary)
summary_result <- to_plot %>%
  group_by(V4, annotation_label) %>%
  summarize(total_count = sum(count))
summary_with_percentage <- summary_result %>%
  group_by(V4) %>%
  mutate(percentage = total_count / sum(total_count) * 100)
to_plot <- as.data.frame(summary_with_percentage)
dictionary <- list("E1"=1, "E2"=2, "E3"=3,
                   "E4"=4, "E5"=5, "E6"=7,
                   "E7"=8, "E8"=6, "E9"=9,
                   "E10"=10,"E11"=11,"E12"=15,
                   "E13"=12,"E14"=13,"E15"=14)
keys <- names(dictionary)
values <- unlist(dictionary)
to_plot$V4 <- values[match(to_plot$V4, keys)]
to_plot$V4 <- paste0("E",to_plot$V4)
to_plot$V4 <- factor(to_plot$V4,levels=rev(paste0("E",1:state_num)))
to_plot$annotation_label <- factor(to_plot$annotation_label,levels = rev(c("Promoter (<=1kb)","Promoter (1-2kb)","Promoter (2-3kb)","Exon","Intron","Downstream (<=300bp)","Distal Intergenic","3' UTR","5' UTR")))
color <- setNames(c("#b83b5e","#f08a5d","#f9ed69","#11999e","#30e3ca","#112d4e","#3f72af","#cca8e9","#c3bef0"),c("Promoter (<=1kb)","Promoter (1-2kb)",
                        "Promoter (2-3kb)","Exon","Intron",
                        "Downstream (<=300bp)","Distal Intergenic",
                        "3' UTR","5' UTR"))
p <- ggplot(to_plot, aes(x = percentage, y = V4, fill = annotation_label)) +  
  geom_bar(stat = 'identity',color="black") +   
  theme_minimal() +   
  scale_fill_manual(values = color)+
  ylab("Proportion")
ggsave("result/Sup_figures/chromHMM_15_state_annotation_old.pdf",p,width = 6,height = 8)
