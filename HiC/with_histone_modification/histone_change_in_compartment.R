rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(data.table)
options(scipen = 999)  
tissue <- "lung"
resolution <- "50000"
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 
antibody <- "H3K27me3"
histone_change_in_compartment <- function(tissue,antibody,resolution){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  histone <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff.csv"))
  histone$Start <- histone$Start+1
  compartment <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_",resolution,".csv"))
  compartment$start <- compartment$start+1
  
  histone_increase <- data.table(histone[which(histone$Significant=="Up"),c("Chr","Start","End")])
  histone_decrease <- data.table(histone[which(histone$Significant=="Down"),c("Chr","Start","End")])
  compartment <- data.table(compartment)
  setDT(compartment)
  setkey(compartment, chr, start, end) 
  
  if(nrow(histone_increase)>0){
    setDT(histone_increase)
    setkey(histone_increase, Chr, Start, End) 
    increase_overlaps <- foverlaps(compartment, histone_increase, type = "any", nomatch = 0L)
    increase_overlaps_compartment_summary <- as.data.frame(table(increase_overlaps$condition))
    colnames(increase_overlaps_compartment_summary)[2] <- "increase"
    increase_overlaps_compartment_summary$increase_percent <- increase_overlaps_compartment_summary$increase/sum(increase_overlaps_compartment_summary$increase)*100
    }
  if(nrow(histone_increase)>0){
    setDT(histone_decrease)
    setkey(histone_decrease, Chr, Start, End) 
    decrease_overlaps <- foverlaps(compartment, histone_decrease, type = "any", nomatch = 0L) 
    decrease_overlaps_compartment_summary <- as.data.frame(table(decrease_overlaps$condition))
    colnames(decrease_overlaps_compartment_summary)[2] <- "decrease"
    decrease_overlaps_compartment_summary$decrease_percent <- decrease_overlaps_compartment_summary$decrease/sum(decrease_overlaps_compartment_summary$decrease)*100
  }
  
  to_plot <- merge(increase_overlaps_compartment_summary[,c(1,3)],decrease_overlaps_compartment_summary[,c(1,3)],by="Var1")
  color <- read.table("data/samples/7_distinct_color.txt")
  color <- setNames(color$V1,c("A-A","B-B","A-B","B-A"))
  to_plot <- reshape2::melt(to_plot)
  to_plot$variable <- as.character(to_plot$variable)
  to_plot$variable[which(to_plot$variable=="increase_percent")] <- "Increase"
  to_plot$variable[which(to_plot$variable=="decrease_percent")] <- "Decrease"
  to_plot$variable <- factor(to_plot$variable,levels = c("Increase","Decrease"))
  ggplot(to_plot, aes(x = variable, y = value, fill = Var1)) +  
    geom_bar(stat = 'identity',colour = "white") +   
    theme_minimal() +   
    scale_fill_manual(values = color) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue)," ",antibody))
  }

antibody <- "H3K27ac"
histone_in_compartment_change <- function(tissue,resolution){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue_label_change(tissue) & search_table$antibody==antibody),]
  search_table$age[which(search_table$age=="3m")] <- "young"
  search_table$age[which(search_table$age=="24m")] <- "old"
  search_table$label <- paste0(search_table$sample_name,".",search_table$age,".",search_table$mouse_ID)
  histone <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins_diff.csv"))
  histone <- histone[,c("Chr","Start","End",search_table$label),]
  histone$Start <- histone$Start +1
  
  compartment <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_",resolution,".csv"))
  compartment$start <- compartment$start+1
  setDT(compartment)
  setkey(compartment, chr, start, end) 
  
  setDT(histone)
  setkey(histone, Chr, Start, End) 

  overlaps <- foverlaps(compartment, histone, type = "any", nomatch = 0L) 
  overlaps$position <- paste0(overlaps$chr,"-",overlaps$Start,"-",overlaps$End)
  
  B2A <- as.data.frame(overlaps[which(overlaps$condition=="B-A")])
  B2A <- B2A[,c("position",search_table$label)]
  A2B <- as.data.frame(overlaps[which(overlaps$condition=="A-B")])
  A2B <- A2B[,c("position",search_table$label)]
  
  B2A <- reshape2::melt(B2A)
  A2B <- reshape2::melt(A2B)
  B2A$condition <- "B-A"
  A2B$condition <- "A-B"
  to_plot <- rbind(A2B,B2A)
  colnames(to_plot)[2] <- "label"
  to_plot <- merge(to_plot,search_table,by="label")
  to_plot$age <- factor(to_plot$age, levels=c("young","old"))
  search_table$age <- factor(search_table$age,c("young","old"))
  search_table <- search_table[order(search_table$age),]
  to_plot$label <- factor(to_plot$label,levels=search_table$label)
  ggplot(to_plot[which(to_plot$condition=="A-B"),], aes(x = label, y = value,fill=age)) +  
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot() +  
    theme_minimal() +
    theme(text = element_text(size = 20), axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(title = paste0(tissue_label_change(tissue),"\n",antibody," in A-B regions"), x = NULL, y = "CPM")
  
  ggplot(to_plot[which(to_plot$condition=="B-A"),], aes(x = label, y = value,fill=age)) +  
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot() +  
    theme_minimal() +
    theme(text = element_text(size = 20), axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(title = paste0(tissue_label_change(tissue),"\n",antibody," in B-A regions"), x = NULL, y = "CPM")
  
  to_plot$logCPM <- log2(to_plot$value)
  ggplot(to_plot[which(to_plot$condition=="A-B"),], aes(x = label, y = logCPM,fill=age)) +  
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot() +  
    theme_minimal() +
    theme(text = element_text(size = 20), axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(title = paste0(tissue_label_change(tissue),"\n",antibody," in A-B regions"), x = NULL, y = "log2(CPM)")
  
  ggplot(to_plot[which(to_plot$condition=="B-A"),], aes(x = label, y = logCPM,fill=age)) +  
    scale_fill_brewer(palette = "Pastel1") +
    geom_boxplot() +  
    theme_minimal() +
    theme(text = element_text(size = 20), axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(title = paste0(tissue_label_change(tissue),"\n",antibody," in B-A regions"), x = NULL, y = "log2(CPM)")
  
  }

H3K9me3_H3K27me3_reverse_in_compartment <- function(tissue,resolution){
  histone <- read.table(paste0("data/samples/",tissue,"/H3K27me3_H3K9me3_intersect/10kb_all_significant_second_quadrant.bed"))
  colnames(histone)[1:3] <- c("Chr","Start","End")
  histone$Start <- histone$Start+1
  
  compartment <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_",resolution,".csv"))
  compartment$start <- compartment$start+1
  setDT(compartment)
  setkey(compartment, chr, start, end) 
  
  setDT(histone)
  setkey(histone, Chr, Start, End) 
  overlaps <- foverlaps(compartment, histone, type = "any", nomatch = 0L) 
  to_plot <- as.data.frame(table(overlaps$condition))
  
  color <- read.table("data/samples/7_distinct_color.txt")
  color <- setNames(color$V1,c("A-A","B-B","A-B","B-A"))
  
  to_plot$percent <- to_plot$Freq/sum(to_plot$Freq) *100
  to_plot$condition <- "H3K9me3 decrease H3K27me3 increase"
  ggplot(to_plot, aes(x = condition, y = percent, fill = Var1)) +  
    geom_bar(stat = 'identity',colour = "white") +   
    theme_minimal() +   
    scale_fill_manual(values = color) +
    theme(axis.title.x = element_blank(),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue)))
}
