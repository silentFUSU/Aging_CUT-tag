rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(limma)
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
tissue <- "kidney"
resolution <- "10000"
insulation_redundant_TAD <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  boundary_list <- data.frame()
  for(sample in search_table$sample_name){
    boundary <- read.table(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",sample,"/",sample,"_",resolution,"_dense.is500001.ids200001.insulation.boundaries"))
    boundary <- boundary[,c(7:8)]
    boundary <- boundary %>%  
      separate(V7, into = c("bin", "org", "pos"), sep = "\\|") %>%  
      separate(pos, into = c("chr", "start_end"), sep = ":") %>%  
      separate(start_end, into = c("start", "end"), sep = "-")  
    boundary <- boundary[,-c(1:2)]
    colnames(boundary)[4] <- "boundaryScore"
    boundary$sample <- sample
    boundary$chr <- factor(boundary$chr, levels = paste0("chr",c(1:19,"X","Y")))
    boundary$start <- as.numeric(boundary$start)
    boundary$end <- as.numeric(boundary$end)
    boundary_list <- rbind(boundary_list,boundary)
  }
  sorted_boundary_list <- boundary_list %>%
    arrange(chr, start, boundaryScore)  
  sorted_boundary_list$choose_or_not <- "No"
  i=1
  redundant_boundary <- data.frame()
  while(i < nrow(boundary_list)){
    j=i+1
    while(sorted_boundary_list[i,"chr"]==sorted_boundary_list[j,"chr"] &  sorted_boundary_list[j,"start"]-sorted_boundary_list[i,"start"] <= 50000 & j <=nrow(boundary_list) ){
        if(sorted_boundary_list[i,"boundaryScore"] < sorted_boundary_list[j,"boundaryScore"]){
          i = j
          j+1
        }else{
          j=j+1
        }
    }
    sorted_boundary_list[i,"choose_or_not"] <- "Yes"
    redundant_boundary <- rbind(redundant_boundary,sorted_boundary_list[i,])
    i=j
  }
  redundant_boundary <- redundant_boundary[,-ncol(redundant_boundary)]
  write.csv(redundant_boundary,paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD.csv"),row.names = F)
  redundant_boundary$label <- paste0(redundant_boundary$chr,":",redundant_boundary$start,"-",redundant_boundary$end)
  redundant_boundary_list <- paste0(redundant_boundary$chr,":",redundant_boundary$start,"-",redundant_boundary$end)
  bedpe <- data.frame()
  for(row_num in c(1:(nrow(redundant_boundary)-1))){
    t_bedpe <- data.frame(chr=redundant_boundary[row_num,1],
                        start = (redundant_boundary[row_num,2]+redundant_boundary[row_num,3])/2, 
                        end = (redundant_boundary[row_num+1,2]+redundant_boundary[row_num+1,3])/2)
    bedpe <- rbind(bedpe,t_bedpe)  
  }
      

  bedpe <- data.frame(chr1=bedpe$chr, x1=bedpe$start, x2=bedpe$end, 
                      chr2=bedpe$chr, y1=bedpe$start, y2=bedpe$end)
  write.table(bedpe,paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/all_samples_",resolution,"_redundant_tads.bedpe"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  }
tissues <- c("kidney","brain","CB","liver","colon")
tissues <- c("stomach","heart","bonemarrow")
for(tissue in tissues){
  insulation_redundant_TAD(tissue,resolution)
}


boundary_strength_caculate <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  boundary_strength <- list()
  for(sample in search_table$sample_name){
    boundary <- read.table(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",sample,"/",sample,"_",resolution,"_dense.is500001.ids200001.insulation"))  
    boundary <- boundary[,c(1,9)]
    boundary <- boundary %>%  
      separate(V1, into = c("bin", "org", "pos"), sep = "\\|") %>%  
      separate(pos, into = c("chr", "start_end"), sep = ":") %>%  
      separate(start_end, into = c("start", "end"), sep = "-")  
    boundary <- boundary[,-c(1:2)]
    colnames(boundary)[4] <- "delta"
    boundary$label <- paste0(boundary$chr,":",boundary$start,"-",boundary$end)
    boundary_strength[[sample]] <- data.frame()
    for(label in redundant_boundary_list){
      site <- which(boundary$label == label) 
      i <- site-1
      left_delta <- boundary[i,"delta"]
      while(boundary[i,"delta"] < boundary[i-1,"delta"] & boundary[i,"chr"] == boundary[i-1,"chr"] & i-1 > 0){
        i<-i-1
        left_delta <- boundary[i,"delta"]
      }
      
      i <- site+1
      right_delta <- boundary[i,"delta"]
      while(boundary[i,"delta"] > boundary[i+1,"delta"] & boundary[i,"chr"] == boundary[i+1,"chr"] & i+1 <= nrow(boundary)){
        i<-i+1
        right_delta <- boundary[i,"delta"]
      }
      
      if(left_delta >=0 & right_delta >= 0){
        delat <- abs(right_delta - left_delta)
      }else if(left_delta >=0 & right_delta < 0){
        delta <- left_delta - right_delta
      }else if(left_delta <0 & right_delta >=0){
        delta <- right_delta - left_delta
      }else{
        delta <- abs(right_delta - left_delta)
      }
      boundary_strength[[sample]] <- rbind(boundary_strength[[sample]], data.frame(label=label,delta=delta))
    }
    colnames(boundary_strength[[sample]])[2] <- sample
  }
  redundant_boundary_table <- Reduce(function(x, y) merge(x, y, by = "label"), boundary_strength)  
  rownames(redundant_boundary_table) <- redundant_boundary_table[,1]
  redundant_boundary_table <- redundant_boundary_table[,-1]
  design <- cbind(Intercept=1,Group=c(0,0,1,1))
  fit <- lmFit(redundant_boundary_table,design)
  fit <- eBayes(fit)
  out <- topTable(fit,coef=2,number = nrow(redundant_boundary_table))
  out$Significant <- ifelse(out$adj.P.Val < 0.05 & abs(out$logFC) >= 0, 
                            ifelse(out$logFC > 0, "Up", "Down"), "Stable")
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  ggplot(
    out, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," boundary strength"))+
    annotate("text", x = min(out$logFC), y = max(-log10(out$adj.P.Val)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$logFC), y = max(-log10(out$adj.P.Val)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  
  } 