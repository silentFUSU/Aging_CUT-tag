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
    combined_matrix <- combined_matrix[which(abs(combined_matrix$region2 - combined_matrix$region1) > 1000000),]
    matrix_list[[sample]] <- combined_matrix
  }
  numCores <- 5
  # register(MulticoreParam(workers = numCores), default = TRUE) 
  hicexp <- make_hicexp(matrix_list[[search_table$sample_name[which(search_table$age=="3M")][1]]],matrix_list[[search_table$sample_name[which(search_table$age=="3M")][2]]],
                        matrix_list[[search_table$sample_name[which(search_table$age=="24M")][1]]],matrix_list[[search_table$sample_name[which(search_table$age=="24M")][2]]],
                        groups =c(0,0,1,1),
                        zero.p = 0.8, A.min = 5, filter = TRUE)
  
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
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/HiCcompare_",resolution,"_VolcanoPlot_filter_distance.png"),p,width = 5,height = 5, type="cairo")
  saveRDS(hicexp,paste0("data/samples/HiC/",tissue,"/HiCcompare_input_",resolution,"_filter_distance.rds"))
  write.csv(out, paste0("data/samples/HiC/",tissue,"/HiCcompare_output_",resolution,"_filter_distance.csv"),row.names = F)
  
}