rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(diffHic)
library(ggplot2)
library(stringr)
library(edgeR)
library(data.table)
library(GenomicRanges) 
library(rtracklayer)  
library(Matrix)
library(csaw)
library(multiHiCcompare)
library(HiCcompare)
library(BiocParallel)
library(dplyr) 
library(tidyverse)
tissue <- "brain"
resolution <- "1000000"
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
HiCcompare_input_convert <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  matrix_list <- list()
  valid_chr <- paste0("chr",c(1:19,"X","Y"))
  for(sample in search_table$sample_name){
    mat <- fread(paste0("data/samples/HiC/",tissue,"/raw_matrix/",sample,"_",resolution,".matrix"))
    mat <- as.data.frame(mat)
    bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",sample,"_",resolution,"_abs.bed"))
    matrix <- hicpro2bedpe(mat, bed)
    matrix <- matrix$cis[names(matrix$cis) %in% valid_chr]  
    selected_matrix <- lapply(matrix, function(df) {  
      df[, c("chr1", "start1", "start2", "IF"), drop = FALSE]
    })  
    combined_matrix <- Reduce(function(x, y) rbind(x, y), selected_matrix)  
    colnames(combined_matrix) <- c("chr","region1","region2","IF")
    matrix_list[[sample]] <- combined_matrix
  }
  numCores <- 5
  register(MulticoreParam(workers = numCores), default = TRUE) 
  hicexp <- make_hicexp(matrix_list[[search_table$sample_name[which(search_table$age=="3M")][1]]],matrix_list[[search_table$sample_name[which(search_table$age=="3M")][2]]],
                        matrix_list[[search_table$sample_name[which(search_table$age=="24M")][1]]],matrix_list[[search_table$sample_name[which(search_table$age=="24M")][2]]],
                        groups =c(0,0,1,1),
                        zero.p = 0.5, A.min = 10, filter = TRUE)
  
  hicexp <- fastlo(hicexp, verbose = T, parallel = FALSE)
  d <- model.matrix(~factor(meta(hicexp)$group))
  hicexp <- hic_glm(hicexp, design = d, coef = 2, method = "QLFTest", p.method = "fdr", parallel = FALSE)
  compairson <- as.data.frame(hicexp@comparison)
  hic_table <- as.data.frame(hicexp@hic_table)
  out <- cbind(hic_table,compairson[,c(5:9)])
  out$Significant <- ifelse(out$p.adj < 0.05 & abs(out$logFC) >= 0, 
                            ifelse(out$logFC > 0, "Up", "Down"), "Stable")
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(
    out, aes(x = logFC, y = -log10(p.adj))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (fdr)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)))+
    annotate("text", x = min(out$logFC), y = max(-log10(out$p.adj)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$logFC), y = max(-log10(out$p.adj)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  
  dir.create(paste0("result/HiC/",tissue),showWarnings = F)
  dir.create(paste0("result/HiC/",tissue,"/differential_analysis/"),showWarnings = F)
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/HiCcompare_",resolution,"_VolcanoPlot.png"),p,width = 5,height = 5, type="cairo")
  dir.create(paste0("data/samples/HiC/",tissue,"/differential_analysis/"),showWarnings = F,recursive = T)
  saveRDS(hicexp,paste0("data/samples/HiC/",tissue,"/differential_analysis/HiCcompare_input_",resolution,".rds"))
  write.csv(out, paste0("data/samples/HiC/",tissue,"/differential_analysis/HiCcompare_output_",resolution,".csv"),row.names = F)
}

HiCcompare_compartment_annotation_larger_resolution <- function(tissue,resolution){
  out <- fread(paste0("data/samples/HiC/",tissue,"/differential_analysis/HiCcompare_output_",resolution,".csv"), sep = ",")
  out <- as.data.frame(out)
  out$chr <- paste0("chr", out$chr)  
  out$chr[which(out$chr=="chr23")] <- "chrX"
  out$chr[which(out$chr=="chr24")] <- "chrY"
  colnames(out)[2:3] <- c("region1_start","region2_start")
  out$region1_end <- out$region1_start+as.numeric(resolution)
  out$region2_end <- out$region2_start+as.numeric(resolution)
  bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/DYQ135_",resolution,"_abs.bed"))
  valid_chr <- paste0("chr",c(1:19,"X","Y"))
  bed <- bed[which(bed$V1 %in% valid_chr),]
  bed$V2 <- bed$V2 + 1
  bed <- as.data.table(bed)
  setDT(bed)
  setkey(bed, V1, V2, V3) 
  
  compartment_list <- list()
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_50000.PC1.txt"))
    df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    df <- df[,c(2:4,6)]
    df$compartment <- ifelse(df[, 4] > 0, "A", "B")
    df$label <- paste0(df$V2,"-",df$V3,"-",df$V4)
    compartment_list[[i]] <- df
    names(compartment_list)[i]<-sample
  }
  selected_compartment_list <- lapply(compartment_list, function(df) {  
    df[, c("label","compartment"), drop = FALSE]
  }) 
  compartment <- Reduce(function(x, y) merge(x, y, by = "label"), selected_compartment_list)
  rows_to_keep <- apply(compartment[, -1], 1, function(row) all(row == row[1])) 
  compartment <- compartment[rows_to_keep,1:2]
  colnames(compartment)[2] <- "compartment"
  compartment <- separate(compartment, col = label, into = c("chr", "start", "end"), sep = "-")  
  compartment$start <- as.numeric(compartment$start)
  compartment$end <- as.numeric(compartment$end)
  compartment$chr <- factor(compartment$chr,levels=valid_chr)
  comparment <- compartment[order(compartment$chr,compartment$start),]
  compartment$start <- compartment$start + 1
  comparment <- as.data.table(compartment)
  setDT(compartment)
  setkey(compartment, chr, start, end) 
  
  overlaps <- foverlaps(bed,compartment, type = "any", nomatch = 0L)
  overlaps$label <- paste(overlaps$V1,overlaps$V2,overlaps$V3,sep = "-")
  group_counts <- overlaps %>%  
    group_by(label, compartment) %>%  
    summarise(count = n(), .groups = "drop")  
  total_counts <- overlaps %>%  
    group_by(label) %>%  
    summarise(total = n(), .groups = "drop")  
  result <- group_counts %>%  
    inner_join(total_counts, by = "label") %>%  
    mutate(percentage = count / total * 100) 
  unique_result <- result %>%  
    group_by(label) %>%  
    slice_max(percentage) %>%  
    ungroup()  
  processed_result <- unique_result %>%  
    group_by(label) %>%  
    mutate(equal_percentage = if_else(n() == 2 & all(percentage == percentage[1]), TRUE, FALSE)) %>%  
    mutate(compartment = if_else(equal_percentage, "weak alternation region", compartment)) %>%  
    filter(!(equal_percentage & compartment != "weak alternation region")) %>%  
    select(-equal_percentage) %>%  
    slice_max(percentage, with_ties = FALSE) %>%  
    ungroup()  
  
  out$region1_start <- out$region1_start +1 
  out$region2_start <- out$region2_start +1
  out$region1_label <- paste(out$chr,out$region1_start,out$region1_end,sep="-")
  out$region2_label <- paste(out$chr,out$region2_start,out$region2_end,sep="-")
  out <- merge(out,processed_result[,c("label","compartment")],  by.x = "region1_label", by.y = "label", all.x=T)
  out <- out %>%  
    rename(region1_compartment = compartment)  
  out$region1_compartment[is.na(out$region1_compartment)] <- "Unknown"  
  
  out <- merge(out,processed_result[,c("label","compartment")],  by.x = "region2_label", by.y = "label", all.x=T)
  out <- out %>%  
    rename(region2_compartment = compartment)  
  out$region2_compartment[is.na(out$region2_compartment)] <- "Unknown"  
  out$condition <- paste(out$region1_compartment,out$region2_compartment,sep='-')
  out$condition[which(out$condition=="B-A")] <- "A-B"
  out <- out[which(out$condition %in% c("A-A","A-B","B-B")),]
  
  conditions <- c("A-A","A-B","B-B")
  condition_change_percent <- data.frame()
  for(condition in conditions){
    t_out <- out[which(out$condition==condition & out$Significant!="Stable"),]  
    t_condition_change_percent <- as.data.frame(table(t_out$Significant))
    sum = sum(t_condition_change_percent$Freq)
    t_condition_change_percent$Freq <- t_condition_change_percent$Freq/sum*100
    t_condition_change_percent <- data.frame(condition=condition,
                                             Up=t_condition_change_percent$Freq[which(t_condition_change_percent$Var1=="Up")],
                                             Down=t_condition_change_percent$Freq[which(t_condition_change_percent$Var1=="Down")])
    condition_change_percent <- rbind(condition_change_percent,t_condition_change_percent)
  }
  to_plot <- reshape2::melt(condition_change_percent)
  ggplot(to_plot, aes(x = condition, y = value, fill = variable)) +  
    geom_bar(stat = 'identity',colour = "white") +   
    theme_minimal() +   
    scale_fill_brewer(palette = "Pastel1") +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20),legend.title = element_blank()) +
    ylab("Proportion")+
    ggtitle(paste0(tissue_label_change(tissue)),"Compartment Interaction")
  }


HiCcompare_compartment_annotation_same_resolution <- function(tissue,resolution){
  out <- fread(paste0("data/samples/HiC/",tissue,"/differential_analysis/HiCcompare_output_",resolution,".csv"), sep = ",")
  out <- as.data.frame(out)
  out$chr <- paste0("chr", out$chr)  
  out$chr[which(out$chr=="chr23")] <- "chrX"
  out$chr[which(out$chr=="chr24")] <- "chrY"
  
  compartment_list <- list()
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_50000.PC1.txt"))
    df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X","Y")))),]
    df <- df[,c(2:4,6)]
    df$compartment <- ifelse(df[, 4] > 0, "A", "B")
    compartment_list[[i]] <- df
    names(compartment_list)[i]<-sample
  }
  out_annotation_list <- list()
  for(i in c(1:length(search_table$sample_name))){
    sample <- search_table$sample_name[i]
    out_label <- out[,c(1:3)]    
    out_label$region1_compartment <- "A"
    out_label$region2_compartment <- "A"
    out_label$region1_label <- paste0(out_label$chr,"-",out_label$region1)
    out_label$region2_label <- paste0(out_label$chr,"-",out_label$region2)
    compartment_list[[sample]]$label <- paste0(compartment_list[[sample]]$V2,"-",compartment_list[[sample]]$V3)
    out_label <- out_label[which(out_label$region1_label %in% compartment_list[[sample]]$label),]
    out_label <- out_label[which(out_label$region2_label %in% compartment_list[[sample]]$label),]
    out_label$region1_compartment[which(out_label$region1_label %in% compartment_list[[sample]]$label[which(compartment_list[[sample]]$compartment=="B")])] <- "B"
    out_label$region2_compartment[which(out_label$region2_label %in% compartment_list[[sample]]$label[which(compartment_list[[sample]]$compartment=="B")])] <- "B"
    out_label$condition <- paste0(out_label$region1_compartment,"-",out_label$region2_compartment)
    out_label <- out_label[,c("region1_label","region2_label","condition")]
    out_label$contact <- paste0(out_label$region1_label,"_",out_label$region2_label)
    out_annotation_list[[sample]] <- out_label
  }
  selected_out_annotation_list <- lapply(out_annotation_list, function(df) {  
    df[, c("contact","condition"), drop = FALSE]
  })  
  # out_annotation <- selected_out_annotation_list[[1]]
  out_annotation <- Reduce(function(x, y) merge(x, y, by = "contact"), selected_out_annotation_list)
  rows_to_keep <- apply(out_annotation[, -1], 1, function(row) all(row == row[1]))  
  out_annotation <- out_annotation[rows_to_keep,]
  out_annotation <- out_annotation[,c(1:2)]
  colnames(out_annotation) <- c("contact","condition")
  out_annotation$condition[which(out_annotation$condition=="B-A")] <- "A-B"
  out$condition <- "Unknown or Mismatch"
  out$contact <- paste0(out$chr,"-",out$region1,"_",out$chr,"-",out$region2)
  out$condition[which(out$contact %in% out_annotation$contact[which(out_annotation$condition=="A-A")])] <- "A-A"
  out$condition[which(out$contact %in% out_annotation$contact[which(out_annotation$condition=="A-B")])] <- "A-B"
  out$condition[which(out$contact %in% out_annotation$contact[which(out_annotation$condition=="B-B")])] <- "B-B"
  colour <- setNames(c("grey","#00adb5","#ff2e63","#f9ed69"),c("Unknown or Mismatch","A-A","A-B","B-B"))
  p <- ggplot(
    out, aes(x = logFC, y = -log10(p.adj))) +
    geom_point(aes(color = condition), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (fdr)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(paste0(tissue_label_change(tissue)))
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/HiCcompare_",resolution,"_VolcanoPlot_contact_condition.png"),p,width = 7,height = 5, type="cairo")
 
   p1<-ggplot() +
    geom_point(out, mapping=aes(x = logFC, y = -log10(p.adj)), color = "grey",size = 2) +
    geom_point(data=out[which(out$condition=="A-B"),], mapping=aes(x = logFC, y = -log10(p.adj)),color = "#ff2e63",alpha=0.7,size = 2)+
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (fdr)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(paste0(tissue_label_change(tissue)," A-B")) +
    annotate("text", x = min(out$logFC), y = max(-log10(out$p.adj)), label = nrow(out[which(out$condition=="A-B" & out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$logFC), y = max(-log10(out$p.adj)), label = nrow(out[which(out$condition=="A-B" & out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
   ggsave(paste0("result/HiC/",tissue,"/differential_analysis/HiCcompare_",resolution,"_VolcanoPlot_contact_condition_A-B.png"),p1,width = 5,height = 5, type="cairo")
  
   p2<-ggplot() +
     geom_point(out, mapping=aes(x = logFC, y = -log10(p.adj)), color = "grey",size = 2) +
     geom_point(data=out[which(out$condition=="A-A"),], mapping=aes(x = logFC, y = -log10(p.adj)),color = "#00adb5",alpha=0.7,size = 2)+
     geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
     geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
     labs(x="log2(fold change)",
          y="-log10 (fdr)") +
     theme_bw()+
     theme(text = element_text(size = 20))+
     ggtitle(paste0(tissue_label_change(tissue)," A-A")) +
     annotate("text", x = min(out$logFC), y = max(-log10(out$p.adj)), label = nrow(out[which(out$condition=="A-A" & out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
     annotate("text", x = max(out$logFC), y = max(-log10(out$p.adj)), label = nrow(out[which(out$condition=="A-A" & out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
   ggsave(paste0("result/HiC/",tissue,"/differential_analysis/HiCcompare_",resolution,"_VolcanoPlot_contact_condition_A-A.png"),p2,width = 5,height = 5, type="cairo")
   
   p3<-ggplot() +
     geom_point(out, mapping=aes(x = logFC, y = -log10(p.adj)), color = "grey",size = 2) +
     geom_point(data=out[which(out$condition=="B-B"),], mapping=aes(x = logFC, y = -log10(p.adj)),color = "#f9ed69",alpha=0.7,size = 2)+
     geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
     geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
     labs(x="log2(fold change)",
          y="-log10 (fdr)") +
     theme_bw()+
     theme(text = element_text(size = 20))+
     ggtitle(paste0(tissue_label_change(tissue)," B-B")) +
     annotate("text", x = min(out$logFC), y = max(-log10(out$p.adj)), label = nrow(out[which(out$condition=="B-B" & out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
     annotate("text", x = max(out$logFC), y = max(-log10(out$p.adj)), label = nrow(out[which(out$condition=="B-B" & out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
   ggsave(paste0("result/HiC/",tissue,"/differential_analysis/HiCcompare_",resolution,"_VolcanoPlot_contact_condition_B-B.png"),p3,width = 5,height = 5, type="cairo")
   
   conditions <- c("A-A","A-B","B-B")
   condition_change_percent <- data.frame()
   for(condition in conditions){
      t_out <- out[which(out$condition==condition & out$Significant!="Stable"),]  
      t_condition_change_percent <- as.data.frame(table(t_out$Significant))
      sum = sum(t_condition_change_percent$Freq)
      t_condition_change_percent$Freq <- t_condition_change_percent$Freq/sum*100
      t_condition_change_percent <- data.frame(condition=condition,
                                               Up=t_condition_change_percent$Freq[which(t_condition_change_percent$Var1=="Up")],
                                               Down=t_condition_change_percent$Freq[which(t_condition_change_percent$Var1=="Down")])
      condition_change_percent <- rbind(condition_change_percent,t_condition_change_percent)
   }
   to_plot <- reshape2::melt(condition_change_percent)
   ggplot(to_plot, aes(x = condition, y = value, fill = variable)) +  
     geom_bar(stat = 'identity',colour = "white") +   
     theme_minimal() +   
     scale_fill_brewer(palette = "Pastel1") +
     theme(axis.title.x = element_blank(), 
           axis.text.x = element_text(angle = 45, hjust = 1),
           text = element_text(size = 20),legend.title = element_blank()) +
     ylab("Proportion")+
     ggtitle(paste0(tissue_label_change(tissue)),"Compartment Interaction")
   return(out)
}

out_filter_plot <- function(out){
  out_filter <- out[which(abs(out$region2 - out$region1) > 1000000),]
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(
    out_filter, aes(x = logFC, y = -log10(p.adj))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (fdr)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)))+
    annotate("text", x = min(out_filter$logFC), y = max(-log10(out_filter$p.adj)), label = nrow(out_filter[which(out_filter$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out_filter$logFC), y = max(-log10(out_filter$p.adj)), label = nrow(out_filter[which(out_filter$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/HiCcompare_",resolution,"_VolcanoPlot_filter_distance.png"),p,width = 5,height = 5, type="cairo")
  colour <- setNames(c("grey","#00adb5","#ff2e63","#f9ed69"),c("Unknown or Mismatch","A-A","A-B","B-B"))
  p <- ggplot(
    out_filter, aes(x = logFC, y = -log10(p.adj))) +
    geom_point(aes(color = condition), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (fdr)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(paste0(tissue_label_change(tissue)))
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/HiCcompare_",resolution,"_VolcanoPlot_contact_condition_filter_distance.png"),p,width = 7,height = 5, type="cairo")
  
  p1<-ggplot() +
    geom_point(out_filter, mapping=aes(x = logFC, y = -log10(p.adj)), color = "grey",size = 2) +
    geom_point(data=out_filter[which(out_filter$condition=="A-B"),], mapping=aes(x = logFC, y = -log10(p.adj)),color = "#ff2e63",alpha=0.7,size = 2)+
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (fdr)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(paste0(tissue_label_change(tissue)," A-B")) +
    annotate("text", x = min(out_filter$logFC), y = max(-log10(out_filter$p.adj)), label = nrow(out_filter[which(out_filter$condition=="A-B" & out_filter$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out_filter$logFC), y = max(-log10(out_filter$p.adj)), label = nrow(out_filter[which(out_filter$condition=="A-B" & out_filter$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/HiCcompare_",resolution,"_VolcanoPlot_contact_condition_A-B_filter_distance.png"),p1,width = 5,height = 5, type="cairo")
  
  p2<-ggplot() +
    geom_point(out_filter, mapping=aes(x = logFC, y = -log10(p.adj)), color = "grey",size = 2) +
    geom_point(data=out_filter[which(out_filter$condition=="A-A"),], mapping=aes(x = logFC, y = -log10(p.adj)),color = "#00adb5",alpha=0.7,size = 2)+
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (fdr)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(paste0(tissue_label_change(tissue)," A-A")) +
    annotate("text", x = min(out_filter$logFC), y = max(-log10(out_filter$p.adj)), label = nrow(out_filter[which(out_filter$condition=="A-A" & out_filter$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out_filter$logFC), y = max(-log10(out_filter$p.adj)), label = nrow(out_filter[which(out_filter$condition=="A-A" & out_filter$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/HiCcompare_",resolution,"_VolcanoPlot_contact_condition_A-A_filter_distance.png"),p2,width = 5,height = 5, type="cairo")
  
  p3<-ggplot() +
    geom_point(out_filter, mapping=aes(x = logFC, y = -log10(p.adj)), color = "grey",size = 2) +
    geom_point(data=out_filter[which(out_filter$condition=="B-B"),], mapping=aes(x = logFC, y = -log10(p.adj)),color = "#f9ed69",alpha=0.7,size = 2)+
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="red",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (fdr)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(paste0(tissue_label_change(tissue)," B-B")) +
    annotate("text", x = min(out_filter$logFC), y = max(-log10(out_filter$p.adj)), label = nrow(out_filter[which(out_filter$condition=="B-B" & out_filter$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out_filter$logFC), y = max(-log10(out_filter$p.adj)), label = nrow(out_filter[which(out_filter$condition=="B-B" & out_filter$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/HiCcompare_",resolution,"_VolcanoPlot_contact_condition_B-B_filter_distance.png"),p3,width = 5,height = 5, type="cairo")
  }
  
convert2bedpe <- function(tissue, resolution){
  out <- fread(paste0("data/samples/HiC/",tissue,"/differential_analysis/HiCcompare_output_",resolution,".csv"),sep=",")
  out$chr <- paste0("chr", out$chr)  
  out$chr[which(out$chr=="chr23")] <- "chrX"
  out$chr[which(out$chr=="chr24")] <- "chrY"
  out_increase <- as.data.frame(out[which(out$Significant=="Up"),c("chr","region1","region2")])
  out_decrease <- as.data.frame(out[which(out$Significant=="Down"),c("chr","region1","region2")])
  out_increase <- data.frame(chr1=out_increase$chr, x1=as.numeric(out_increase$region1), x2=c(as.numeric(out_increase$region1)+as.numeric(resolution)), 
                             chr2=out_increase$chr, y1=as.numeric(out_increase$region2), y2=c(as.numeric(out_increase$region2)+as.numeric(resolution)))
  out_decrease <- data.frame(chr1=out_decrease$chr, x1=as.numeric(out_decrease$region1), x2=c(as.numeric(out_decrease$region1)+as.numeric(resolution)), 
                             chr2=out_decrease$chr, y1=as.numeric(out_decrease$region2), y2=c(as.numeric(out_decrease$region2)+as.numeric(resolution)))
  
  write.table(out_increase,paste0("data/samples/HiC/",tissue,"/differential_analysis/HiCcompare_output_",resolution,"_increase.bedpe"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(out_decrease,paste0("data/samples/HiC/",tissue,"/differential_analysis/HiCcompare_output_",resolution,"_decrease.bedpe"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  
  }
