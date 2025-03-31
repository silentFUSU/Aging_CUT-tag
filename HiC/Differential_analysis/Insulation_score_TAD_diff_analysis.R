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
resolution <- "20000"
tissue <- "kidney"
TAD_diff_analysis <- function(tissue,resolution){
  TAD <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD.csv"))
  TAD <- TAD[,c(1:3)]
  TAD_regions <- data.frame(chr=as.character(),start=as.numeric(),end=as.numeric())
  i=2
  while(i <= nrow(TAD)){
    if(TAD[i,"chr"] == TAD[(i-1),"chr"]){
      t_TAD_regions <- data.frame(chr=TAD$chr[i],start=TAD$start[i-1],end=TAD$start[i]) 
      TAD_regions <- rbind(TAD_regions,t_TAD_regions)
      i <- i+1
    }else{
      i <- i+1
    }
  }
  TAD_regions$start <- TAD_regions$start+1
  TAD_regions <- as.data.table(TAD_regions)
  setDT(TAD_regions)
  setkey(TAD_regions,chr,start,end)
  
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue),]
  
  bed <- read.table(paste0("data/samples/HiC/",tissue,"/raw_matrix/",search_table$sample_name[1],"_",resolution,"_abs.bed"))
  bed <- bed[which(bed$V1 %in% paste0("chr",c(1:19,"X"))),]
  bed$V5 <- "NA"
  count <- 1
  bed[1,"V5"] <- count
  for(i in c(2:nrow(bed))){
    if(bed[i,"V1"] != bed[i-1,"V1"]){
      count <- 1
      bed[i,"V5"] <- count
    }else{
      count <- count+1
      bed[i,"V5"] <- count
    }
  }
  
  bed$V2 <- bed$V2+1
  bed <- as.data.table(bed)
  setDT(bed)
  setkey(bed,V1,V2,V3)
  
  overlaps <- foverlaps(TAD_regions, bed, type = "any", nomatch = 0L)  
  overlaps$TAD_label <- paste(overlaps$chr,overlaps$start,overlaps$end,sep = "-")
  colnames(overlaps)[which(colnames(overlaps)=="V4")] <- "bin_label"
  TAD_counts <- data.frame()
  samples <- search_table$sample_name
  for(sample in samples){
    counts <- fread(paste0("data/samples/HiC/",tissue,"/raw_matrix/",sample,"_",resolution,".matrix"))
    filtered_counts <- counts[V1 %in% overlaps$bin_label]  
    filtered_counts <- filtered_counts[V2 %in% overlaps$bin_label]  
    # filtered_counts[, interaction_label := paste0(V2, "-", V1)]  
    # filtered_counts <- filtered_counts[interaction_label %in% re_sig$interaction_label]
    data <- merge(filtered_counts,overlaps[,c("bin_label","TAD_label")],by.x="V1",by.y="bin_label")
    data <- merge(data,overlaps[,c("bin_label","TAD_label")],by.x="V2",by.y="bin_label")
    data <- data[which(data$TAD_label.x == data$TAD_label.y),]
    data <- data[which(abs(data$V1-data$V2)> floor(50000/as.numeric(resolution)) ),]
    data$TAD_label <- data$TAD_label.x
    sum_by_TAD <- data %>%  
      group_by(TAD_label) %>%  
      summarise(total_V3 = sum(V3))  
    colnames(sum_by_TAD)[2] <- sample
    if(nrow(TAD_counts)==0){
      TAD_counts <- sum_by_TAD
    }else{
      TAD_counts <- merge(TAD_counts,sum_by_TAD,by="TAD_label")
    }
  }
  rownames(TAD_counts) <- TAD_counts$TAD_label
  TAD_counts <- TAD_counts[,-1]
  age <- search_table$age
  y= DGEList(counts=TAD_counts,group=age)
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
    ggtitle(paste0(tissue_label_change(tissue)," TAD"))+
    annotate("text", x = min(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old-young`), y = max(-log10(out$`FDR.old-young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_",resolution,"_TAD_condition_edger_volcano.png"),p,width = 5,height = 6,type="cairo")
  write.csv(out,paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD_diff.csv"))
  # re <-  read.table(paste0("data/samples/HiC/",tissue,"/differential_analysis/",tissue,"_",resolution,".FDR"))
  # re$Significant <- "Stable"
  # re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 < 0)] <- "Down"
  # re$Significant[which((re$V5 < 0.01 & re$V6 < 0.01) & re$V4 > 0)] <- "Up"
  # re_sig <- re[which(re$Significant!="Stable"),]
  # re_sig <- re_sig[which(abs(re_sig$V2 - re_sig$V3)>4),]
  # re_sig$V1 <- paste0("chr",re_sig$V1)
  # re_sig$V1[which(re_sig$V1=="chr20")] <- "chrX"
  # 
  # re_sig$label1 <- paste(re_sig$V1,re_sig$V2,sep = "-")
  # re_sig$label2 <- paste(re_sig$V1,re_sig$V3,sep = "-")
  # overlaps$label <- paste(overlaps$chr,overlaps$V5,sep = "-")
  # 
  # re_sig <- merge(re_sig,overlaps[,c("label","TAD_label")],by.x="label1",by.y="label")
  # colnames(re_sig)[which(colnames(re_sig)=="TAD_label")] <- "TAD_label1"
  # 
  # re_sig <- merge(re_sig,overlaps[,c("label","TAD_label")],by.x="label2",by.y="label")
  # colnames(re_sig)[which(colnames(re_sig)=="TAD_label")] <- "TAD_label2"
  # re_sig <- re_sig[which(re_sig$TAD_label1==re_sig$TAD_label2),]
  # result <- re_sig %>%  
  #   group_by(TAD_label1) %>%  
  #   summarise(  
  #     Up = sum(Significant == "Up"),  
  #     Down = sum(Significant == "Down")  
  #   )  
  # result$condition <- "Stable"
  # result$condition[which(result$Up - result$Down >=5)] <- "Up"
  # result$condition[which(result$Down - result$Up >=5)] <- "Down"
  # color <- setNames(c("red","blue","grey"),c("Up","Down","Stable"))
  # out$TAD_label <- rownames(out)
  # result <- merge(result,out[,c("TAD_label","Significant")],by.x="TAD_label1",by.y="TAD_label")
  # p <- ggplot(result, aes(x = Down, y = Up, color=Significant)) +      
  #   geom_point(size = 1) +  
  #   scale_color_manual(values = color) +
  #   labs(title = paste0(tissue_label_change(tissue)," Old vs. Young"),  
  #        x = "# of decreased interactions in TADs",  
  #        y = "# of increased interactions in TADs",
  #        color = "TAD condition") +        
  #   xlim(0,max(result$Down)+10) +
  #   ylim(0,max(result$Up)+10) +
  #   theme_minimal()+
  #   theme(  
  #     # 修改所有字体的大小  
  #     plot.title = element_text(size = 16),  # 图标题字体大小  
  #     axis.title = element_text(size = 14),  # 轴标题字体大小  
  #     axis.text = element_text(size = 12),   # 坐标刻度标签字体大小  
  #     legend.text = element_text(size = 12), # 图例字体大小  
  #     legend.title = element_text(size = 14) # 图例标题字体大小  
  #   )  
  # ggsave(paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_",resolution,"_TAD_condition_edger.png"),p,width = 5,height = 4,type="cairo")
  # 
}
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus")

for(tissue in tissues){
  TAD_diff_analysis(tissue,resolution)
}
