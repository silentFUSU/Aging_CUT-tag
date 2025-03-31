rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(limma)
library(data.table)
library(edgeR)
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

resolution <- "25000"
tissue <- "lung"
loop_diff_analysis <- function(tissue,resolution){
  loop <- read.table(paste0("data/samples/HiC/",tissue,"/loop/HiCCUPS/merged_",resolution,"_loop.bed"),skip = 1,header = F)
  loop <- loop[,c(1:6)]
  loop <- loop[which(loop$V1 %in% paste0("chr",c(1:19,"X")) & loop$V4 %in% paste0("chr",c(1:19,"X"))),]
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue),]
  bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",search_table$sample_name[1],"_10000_abs.bed"))
  bed <- bed[which(bed$V1 %in% paste0("chr",c(1:19,"X"))),]
  
  loop$V2 <- loop$V2 + 1
  loop$V5 <- loop$V5 + 1
  bed$V2 <- bed$V2 + 1
  
  loop$label1 <- paste(loop$V1,loop$V2,loop$V3,sep = "-")
  loop$label2 <- paste(loop$V4,loop$V5,loop$V6,sep = "-")
  bed$label <- paste(bed$V1,bed$V2,bed$V3,sep = "-")
  
  loop <- as.data.table(loop)
  bed <- as.data.table(bed)
  
  setDT(loop)
  setkey(loop,V1,V2,V3)
  setDT(bed)
  setkey(bed,V1,V2,V3)
  
  overlaps <- foverlaps(loop, bed, type = "any", nomatch = 0L)  
  overlaps <- overlaps[,c("V1","V2","V3","V4","label","V5","V6","label1","label2")]
  overlaps <- as.data.table(overlaps)
  setDT(overlaps)
  setkey(overlaps,V1,V5,V6)
  overlaps <- foverlaps(overlaps,bed, type = "any", nomatch = 0L)
  overlaps <- overlaps[,c("label1","label2","V4","i.V4")]
  colnames(overlaps)[c(3,4)] <- c("bin1","bin2")
  overlaps$label <- apply(overlaps, 1, function(row) {  
    bin1 <- as.numeric(row["bin1"])  
    bin2 <- as.numeric(row["bin2"])  
    if (bin1 > bin2) {  
      paste(bin2, bin1, sep = "-")  
    } else {  
      paste(bin1, bin2, sep = "-")  
    }  
  })  
  samples <- search_table$sample_name
  overlaps$loop_label <- paste(overlaps$label1,overlaps$label2,sep = "-")
  Loop_count <- data.frame()
  for(sample in samples){
    counts <- fread(paste0("data/samples/HiC/",tissue,"/raw_matrix/",sample,"_10000.matrix"))
    filtered_counts <- counts[V1 %in% overlaps$bin1 | V1 %in% overlaps$bin2]  
    filtered_counts <- filtered_counts[V2 %in% overlaps$bin1 | V2 %in% overlaps$bin2]
    filtered_counts$label <- paste(filtered_counts$V1,filtered_counts$V2,sep = "-")
    # filtered_counts[, interaction_label := paste0(V2, "-", V1)]  
    # filtered_counts <- filtered_counts[interaction_label %in% re_sig$interaction_label]
    data <- merge(filtered_counts[,c("label","V3")],overlaps,by="label")
    sum_by_Loop <- data %>%  
      group_by(loop_label) %>%  
      summarise(total_V3 = sum(V3))  
    colnames(sum_by_Loop)[2] <- sample
    if(nrow(Loop_count)==0){
      Loop_count <- sum_by_Loop
    }else{
      Loop_count <- merge(Loop_count,sum_by_Loop,by="loop_label")
    }
  }
  rownames(Loop_count) <- Loop_count$loop_label
  Loop_count <- Loop_count[,-1]
  age <- search_table$age
  y= DGEList(counts=Loop_count,group=age)
  y$samples$year <- age
  y$samples$year <- factor(y$samples$year,c("3M","24M"))
  y <- calcNormFactors(y)
  design <- model.matrix(~year, y$samples)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = 2)
  out = cbind(cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
              "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC.old-young"=lrt$table$logFC)
  out <- as.data.frame(out)
  out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= 0, 
                            ifelse(out$`LogFC.old-young` > 0, "Up", "Down"), "Stable")
  write.csv(out,paste0("data/samples/HiC/",tissue,"/loop/HiCCUPS/diff_interaction_within_",resolution,"_loop.csv"))
  colour <- setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(
    out, aes(x = `LogFC.old-young`, y = -log10(`FDR.old-young`))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (FDR)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," Loop"))+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/",tissue,"_HiCCUP_Loop_",resolution,"_diff_volcano_edger.png"),p,width = 5,height = 6,type="cairo")
}
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus")

for(tissue in tissues){
  loop_diff_analysis(tissue,resolution)
}
