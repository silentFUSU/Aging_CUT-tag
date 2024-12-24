rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggrepel)
library(maditr)
library(ChIPseeker)
library(genomation)
library(GenomeInfoDb)
library("GenomicRanges")
state_num <- 14
state <- "E2"
target_state <- "E6"
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
tissues <-  c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
              "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
color_tissues <- sort(c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","jejunum","kidney","liver",
                        "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT","ileum"))

color <- read.table("data/samples/30_distinct_color.txt")
color <- setNames(color$V1,sapply(color_tissues, tissue_label_change))

specific_state_transfer_in_all_tissues <- function(tissues,state_num,state,target_state){
  dir.create(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/",paste0(state, collapse="-"),"_to_",target_state))
  dir.create(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/",paste0(state, collapse="-"),"_to_",target_state,"/bed"))
  transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))
  transfer_matrix <- transfer_matrix %>%  
    group_by(tissue, young_state) %>% 
    mutate(total_freq = sum(Freq), 
           freq_ratio = Freq / total_freq * 100) %>% 
    ungroup()
  transfer_matrix <- transfer_matrix[which(transfer_matrix$young_state==state & transfer_matrix$old_state==target_state),]
  transfer_matrix <- transfer_matrix[order(transfer_matrix$freq_ratio,decreasing = T),]
  transfer_matrix$tissue <- factor(transfer_matrix$tissue, levels=transfer_matrix$tissue)
  ggplot(transfer_matrix, aes(x = tissue, y = freq_ratio,color=tissue)) +  
    geom_point(size = 3) +  
    scale_color_manual(values = color)+
    labs(title = "点图示例", x = NULL, y = "Percentage(%)") +  
    theme_bw()+
    ylim(0,30)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 18))+
    ggtitle(paste0(state," to ",target_state))
  
  transfer_matrix <- read.csv(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/state_transfer/state_transfer.csv"))
  transfer_matrix <- transfer_matrix %>%  
    group_by(tissue) %>% 
    mutate(  
      total_freq = sum(Freq),  
      freq_ratio = Freq / total_freq * 100  
    ) %>%  
    ungroup()
  transfer_matrix <- transfer_matrix[which(transfer_matrix$young_state==state & transfer_matrix$old_state==target_state),]
  
  transfer_matrix <- transfer_matrix[order(transfer_matrix$freq_ratio,decreasing = T),]
  transfer_matrix$tissue <- factor(transfer_matrix$tissue, levels=transfer_matrix$tissue)
  ggplot(transfer_matrix, aes(x = tissue, y = freq_ratio,color=tissue)) +  
    geom_point(size = 3) +  
    scale_color_manual(values = color)+
    labs(title = "点图示例", x = NULL, y = "Percentage(%)") +  
    theme_bw()+
    ylim(0,10)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 18))+
    ggtitle(paste0(state," to ",target_state))
  
}

same_state_transfer_region_in_most_tissues <- function(tissues, state_num, state, target_state){
  state_transfer <- data.frame()
  for(i in c(1:length(tissues))){
    tissue <- tissues[i]
    young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
    young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
    young1$label <- paste(young1$V1,young1$V2,young1$V3,young1$V4,sep = "-")
    young2$label <- paste(young2$V1,young2$V2,young2$V3,young2$V4,sep = "-")
    young <- young1[which(young1$label %in% intersect(young1$label,young2$label)),]
    
    old1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old1_",state_num,"_segments_1k.bed"),header = F)
    old2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old2_",state_num,"_segments_1k.bed"),header = F)
    old1$label <- paste(old1$V1,old1$V2,old1$V3,old1$V4,sep = "-")
    old2$label <- paste(old2$V1,old2$V2,old2$V3,old2$V4,sep = "-")
    old <- old1[which(old1$label %in% intersect(old1$label,old2$label)),]
    
    young$label <- paste0(young$V1,"-",young$V2,"-",young$V3)
    old$label <- paste0(old$V1,"-",old$V2,"-",old$V3)
    young <- young[which(young$label %in% old$label),]  
    old <- old[which(old$label %in% young$label),]
    
    young <- young[which(young$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
    old <- old[which(old$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
    colnames(young)[4] <- "young_state"
    colnames(old)[4] <- "old_state"
    df <- merge(young[,c("label","young_state")],old[,c("label","old_state")],by="label")
    df <- df[which(df$young_state==state & df$old_state==target_state),]
    df$tissue <- tissue
    state_transfer <- rbind(state_transfer,df)
  }
  state_count <- state_transfer %>%   
    count(label)
  state_tissue <- state_transfer %>%   
    group_by(label) %>%   
    summarise(tissue_content = paste(unique(tissue), collapse = "/"))  
  state_count <- merge(state_count,state_tissue,by="label")
}



histone <- "H3K4me3"
condition <- "Up"
specific_state_transfer_in_tissue <- function(tissue,state_num,state,target_state,histone,condition){
  young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
  young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
  young1$label <- paste(young1$V1,young1$V2,young1$V3,young1$V4,sep = "-")
  young2$label <- paste(young2$V1,young2$V2,young2$V3,young2$V4,sep = "-")
  young <- young1[which(young1$label %in% intersect(young1$label,young2$label)),]
  
  old1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old1_",state_num,"_segments_1k.bed"),header = F)
  old2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old2_",state_num,"_segments_1k.bed"),header = F)
  old1$label <- paste(old1$V1,old1$V2,old1$V3,old1$V4,sep = "-")
  old2$label <- paste(old2$V1,old2$V2,old2$V3,old2$V4,sep = "-")
  old <- old1[which(old1$label %in% intersect(old1$label,old2$label)),]
  
  young$label <- paste0(young$V1,":",young$V2,"-",young$V3)
  old$label <- paste0(old$V1,":",old$V2,"-",old$V3)
  young <- young[which(young$label %in% old$label),]  
  old <- old[which(old$label %in% young$label),]
  colnames(young)[4] <- "young_state"
  colnames(old)[4] <- "old_state"
  state_change <- merge(young[,c(4:5)],old[,c(4:5)],by="label")
  state_change <- state_change[which(state_change$young_state==state & state_change$old_state==target_state),]
  if(histone %in% c("H3K27ac","H3K4me3","H3K4me1")){
    bin_size <- "1kb"
  }else{
    bin_size <- "10kb"
  }
  histone_change <- read.csv(paste0("data/samples/",tissue,"/",histone,"/",histone,"_",bin_size,"_bins_diff.csv"))
  histone_change <- histone_change[which(histone_change$Significant==condition),]
  histone_change <- data.frame(label=paste0(histone_change$Chr,":",histone_change$Start,"-",histone_change$End),fdr=histone_change$FDR.old.young,logcpm=histone_change$logCPM)
  state_change <- merge(state_change,histone_change,by="label")
  }

state_num <- 14
state <- "E2"
target_state <- "E6"
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
type1 <- c("liver","bonemarrow","heart","skin","spleen","cecum","colon","lung","brain","Hip","aorta","muscle","stomach")
type2 <- c("ovary","mammarygland","tongue","uterus","thymus","jejunum","testis","iWAT","BAT","kidney","CB","bladder","pancreas")
specific_state_transfer_in_tissue_annotation <- function(tissues,state_num,state,target_state,histone,condition){
  anno_summary <- data.frame()
  for(tissue in tissues){
    young1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young1_",state_num,"_segments_1k.bed"),header = F)
    young2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_young2_",state_num,"_segments_1k.bed"),header = F)
    young1$label <- paste(young1$V1,young1$V2,young1$V3,young1$V4,sep = "-")
    young2$label <- paste(young2$V1,young2$V2,young2$V3,young2$V4,sep = "-")
    young <- young1[which(young1$label %in% intersect(young1$label,young2$label)),]
    
    old1 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old1_",state_num,"_segments_1k.bed"),header = F)
    old2 <- read.delim(paste0("result/all/ChromHMM/all_tissues/",state_num,"_all_tissues/split_1k/",tissue,"_old2_",state_num,"_segments_1k.bed"),header = F)
    old1$label <- paste(old1$V1,old1$V2,old1$V3,old1$V4,sep = "-")
    old2$label <- paste(old2$V1,old2$V2,old2$V3,old2$V4,sep = "-")
    old <- old1[which(old1$label %in% intersect(old1$label,old2$label)),]
    
    young$label <- paste0(young$V1,":",young$V2,"-",young$V3)
    old$label <- paste0(old$V1,":",old$V2,"-",old$V3)
    young <- young[which(young$label %in% old$label),]  
    old <- old[which(old$label %in% young$label),]
    colnames(young)[4] <- "young_state"
    colnames(old)[4] <- "old_state"
    state_change <- merge(young,old[,c(4:5)],by="label")
    state_change <- state_change[which(state_change$young_state==state & state_change$old_state==target_state),]
    state_change <- state_change[which(state_change$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
    state_change_Granges <- GRanges(seqnames = state_change$V1,   
                       ranges = IRanges(start = state_change$V2, end = state_change$V3))
    state_change_anno <- annotatePeak(state_change_Granges, tssRegion=c(-3000, 3000),
                                  TxDb=txdb, annoDb="org.Mm.eg.db")
    state_change_anno_summary <- state_change_anno@annoStat
    state_change_anno_summary$tissue <- tissue_label_change(tissue)
    
    exon_sum <- state_change_anno_summary %>%
      filter(grepl("Exon",Feature)) %>%
      summarise(Frequency = sum(Frequency)) 
    new_exon_row <- data.frame(Feature = "Exon", Frequency = exon_sum$Frequency, tissue=tissue_label_change(tissue))  
    
    intron_sum <- state_change_anno_summary %>%
      filter(grepl("Intron",Feature)) %>%
      summarise(Frequency = sum(Frequency)) 
    new_intron_row <- data.frame(Feature = "Intron", Frequency = intron_sum$Frequency,tissue=tissue_label_change(tissue))  
    
    prom_sum <- state_change_anno_summary %>%
      filter(grepl("Promoter",Feature)) %>%
      summarise(Frequency = sum(Frequency)) 
    Distal_prom <-  data.frame(Feature = "Distal prom", Frequency = (prom_sum$Frequency - state_change_anno_summary$Frequency[which(state_change_anno_summary$Feature == "Promoter (<=1kb)")]), tissue=tissue_label_change(tissue))  
    
    state_change_anno_summary <- rbind(state_change_anno_summary,new_exon_row)
    state_change_anno_summary <- rbind(state_change_anno_summary,new_intron_row)
    state_change_anno_summary <- rbind(state_change_anno_summary,Distal_prom)
    if(nrow(anno_summary)==0){
      anno_summary <- state_change_anno_summary
    }else{
      anno_summary <- rbind(anno_summary,state_change_anno_summary)
    }
  }
  anno_summary$Feature <- as.character(anno_summary$Feature)
  anno_summary$Feature[which(anno_summary$Feature=="Downstream (<=300)")] <- "Downstream"
  anno_summary$Feature[which(anno_summary$Feature=="Distal Intergenic")] <- "Intergenic"
  anno_summary$Feature[which(anno_summary$Feature=="Promoter (<=1kb)")] <- "Promoter"
  anno_to_plot <- anno_summary
  anno_to_plot <- anno_to_plot[which(anno_to_plot$Feature %in% c("Intergenic", "Downstream", "3' UTR", "Intron", "Exon", "5' UTR", "Promoter", "Distal prom")),]
  anno_to_plot$Feature <- factor(anno_to_plot$Feature, levels = c("Intergenic", "Downstream", "3' UTR", "Intron", "Exon", "5' UTR", "Promoter", "Distal prom"))
  color <- setNames(c("#838B8B","#E69900","#BF9915","#99992A","#739940","#4C9955","#26996A","#009980"),c("Intergenic", "Downstream", "3' UTR", "Intron", "Exon", "5' UTR", "Promoter", "Distal prom"))
  ggplot(anno_to_plot, aes(x = tissue, y = Frequency, fill = Feature)) +  
    geom_bar(stat = 'identity',colour = "white") +   
    theme_minimal() +   
    scale_fill_manual(values = color) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Proportion")+
    ggtitle(paste0("Type2 tissues E11 -> E14 Gene location"))
}


state_change_anno_df <- unique(as.data.frame(state_change_anno))
state_change_anno_df_genelist <- bitr(state_change_anno_df$SYMBOL[which(str_detect(state_change_anno_df$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
state_change_anno_df_GO <- enrichGO( state_change_anno_df_genelist$ENTREZID,#GO富集分析
                            OrgDb = GO_database,
                            keyType = "ENTREZID",#设定读取的gene ID类型
                            ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                            pvalueCutoff = 0.05,#设定p值阈值
                            qvalueCutoff = 0.05,#设定q值阈值
                            readable = T)
state_change_anno_df_GO_result <- state_change_anno_df_GO@result
barplot(state_change_anno_df_GO,title = paste0(tissue_label_change(tissue)," E2 -> E6 GO pathway"),label_format = 50)
