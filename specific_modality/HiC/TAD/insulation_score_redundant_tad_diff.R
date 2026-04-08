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
tissue <- "cecum"
resolution <- "20000"

boundary_strength_caculate <- function(tissue,resolution){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  redundant_boundary <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD.csv"))
  redundant_boundary <- redundant_boundary[,c(1:3)]
  colnames(redundant_boundary) <- paste0("redundant_",colnames(redundant_boundary))
  redundant_boundary$label <- paste0(redundant_boundary$redundant_chr,":",redundant_boundary$redundant_start,"-",redundant_boundary$redundant_end)
  boundary_strength <- list()
  for(sample in search_table$sample_name){
    boundary <- read.table(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",sample,"/",sample,"_",resolution,"_dense.is500001.ids200001.insulation.boundaries"))
    boundary <- boundary[,c(7:8)]
    boundary <- boundary %>%  
      separate(V7, into = c("bin", "org", "pos"), sep = "\\|") %>%  
      separate(pos, into = c("chr", "start_end"), sep = ":") %>%  
      separate(start_end, into = c("start", "end"), sep = "-")  
    boundary <- boundary[,-c(1:2)]
    colnames(boundary)[4] <- "boundaryScore"
    boundary$chr <- factor(boundary$chr, levels = paste0("chr",c(1:19,"X","Y")))
    boundary$start <- as.numeric(boundary$start)
    boundary$end <- as.numeric(boundary$end)

    boundary_strength[[sample]] <- data.frame()
    df <- redundant_boundary
    df <- df %>%  
      rowwise() %>%  
      mutate(boundaryScore = {  
        match <- boundary %>%  
          filter(  
            chr == redundant_chr &   
              abs(start - redundant_start) <= 50000
          )  
        if (nrow(match) > 0) {  
          match$boundaryScore[1]
        } else {  
          NA  
        }  
      }) 
    
    boundary <- read.table(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",sample,"/",sample,"_",resolution,"_dense.is500001.ids200001.insulation"))  
    boundary <- boundary[,c(1,9)]
    boundary <- boundary %>%  
      separate(V1, into = c("bin", "org", "pos"), sep = "\\|") %>%  
      separate(pos, into = c("chr", "start_end"), sep = ":") %>%  
      separate(start_end, into = c("start", "end"), sep = "-")  
    boundary <- boundary[,-c(1:2)]
    colnames(boundary)[4] <- "delta"
    boundary$label <- paste0(boundary$chr,":",boundary$start,"-",boundary$end)
    
    for(label in df$label[which(is.na(df$boundaryScore))]  ){
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
        delta <- abs(right_delta - left_delta)
      }else if(left_delta >=0 & right_delta < 0){
        delta <- left_delta - right_delta
      }else if(left_delta <0 & right_delta >=0){
        delta <- right_delta - left_delta
      }else{
        delta <- abs(right_delta - left_delta)
      }
      df[which(df$label==label),"boundaryScore"] <- delta
    }
    boundary_strength[[sample]] <- as.data.frame(df[,c("label","boundaryScore")])
    colnames(boundary_strength[[sample]])[2] <- sample
  }
  redundant_boundary_table <- Reduce(function(x, y) merge(x, y, by = "label"), boundary_strength)  
  rownames(redundant_boundary_table) <- redundant_boundary_table[,1]
  redundant_boundary_table <- redundant_boundary_table[,-1]
  search_table$age[which(search_table$age=="3M")] <- "young"
  search_table$age[which(search_table$age=="24M")] <- "old"
  search_table$sample_name <- factor(search_table$sample_name, levels=colnames(redundant_boundary_table))
  search_table <- search_table[order(search_table$sample_name),]
  search_table$age <- factor(search_table$age,levels=c("young","old"))
  design <- model.matrix(~age, search_table)
  
  fit <- lmFit(redundant_boundary_table,design)
  fit <- eBayes(fit)
  out <- topTable(fit,coef=2,number = nrow(redundant_boundary_table))
  out <- merge(redundant_boundary_table,out,by = "row.names")
  out$Significant <- ifelse(out$adj.P.Val < 0.05 & abs(out$logFC) >= 0, 
                            ifelse(out$logFC > 0, "Up", "Down"), "Stable")
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(
    out, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (FDR)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," boundary strength"))+
    annotate("text", x = min(out$logFC), y = max(-log10(out$adj.P.Val)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$logFC), y = max(-log10(out$adj.P.Val)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  write.csv(out,paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD_diff_boundary.csv"),row.names = F)
  return(p)
} 

tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus")
p_list <- list()
for(tissue in tissues){
  p_list[[tissue]] <- boundary_strength_caculate(tissue,resolution)
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_plot <- plot_a_list(p_list,no_of_rows = 3,no_of_cols = 4)
ggsave(paste0("result/HiC/all_tissues_insulation_score_",resolution,"_boundary_diff.png"),combined_plot,width = 20,height = 15,type="cairo")
