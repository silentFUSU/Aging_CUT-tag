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
tissue <- "ileum"
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
  TAD_regions$length <- TAD_regions$end - TAD_regions$start +1
  TAD_regions <- TAD_regions[which(TAD_regions$length >= 250000),]
  TAD_regions$start <- TAD_regions$start+1
  TAD_regions <- as.data.table(TAD_regions[,c(1:3)])
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
  keep = which(rowSums(cpm(y)>1)>=2)
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
  ggsave(paste0("result/HiC/",tissue,"/differential_analysis/with_histone/",tissue,"_",resolution,"_TAD_condition_edger_volcano_length_larger_250000.png"),p,width = 5,height = 6,type="cairo")
  write.csv(out,paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD_diff_larger_250000.csv"))
}
tissues <- sort(c("brain","CB", "kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip", "mammarygland", "stomach", "thymus","skin","muscle","cecum","ileum"))
tissues <- c("spleen","pancreas")
for(tissue in tissues){
  TAD_diff_analysis(tissue,resolution)
}

p_list <- list()
for(tissue in tissues){
  out <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD_diff_larger_250000.csv"))
  p_list[[tissue]] <- ggplot(
    out, aes(x = `LogFC.old.young`, y = -log10(`FDR.old.young`))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (FDR)") +
    theme_bw()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ggtitle(paste0(tissue_label_change(tissue)," TAD"))+
    annotate("text", x = min(out$`LogFC.old.young`), y = max(-log10(out$`FDR.old.young`)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$`LogFC.old.young`), y = max(-log10(out$`FDR.old.young`)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
combined_plot <- plot_a_list(p_list,no_of_rows = 4,no_of_cols = 4)
ggsave(paste0("result/HiC/all_tissues_insulation_score_",resolution,"_within_TAD_diff_change.png"),combined_plot,width = 18,height = 24,type="cairo")






# summary <- data.frame()
# gap_summary <- data.frame()
# for(tissue in tissues){
#   df <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_",resolution,"_TAD_diff.csv"))
#   t_summary <- data.frame(condition=c("Up","Down"),Freq=c(nrow(df[which(df$Significant=="Up"),]), nrow(df[which(df$Significant=="Down"),])))
#   t_summary$tissue <- tissue_label_change(tissue)
#   t_gap_summary <- data.frame(gap=(nrow(df[which(df$Significant=="Up"),])-nrow(df[which(df$Significant=="Down"),])),tissue=tissue_label_change(tissue))
#   summary <- rbind(summary,t_summary) 
#   gap_summary <- rbind(gap_summary,t_gap_summary)
# }
# gap_summary <- gap_summary[order(gap_summary$gap),]
# summary$tissue <- factor(summary$tissue,levels=gap_summary$tissue)
# 
# 
# summary$Freq[which(summary$condition=="Down")] <- -summary$Freq[which(summary$condition=="Down")]
# to_plot <- summary[-which(summary$tissue %in% c("Colon","Bone Marrow","Liver")),]
# ggplot(to_plot, aes(x = Freq, y = tissue, fill = condition)) +  
#   geom_bar(stat = "identity")+
#   scale_x_continuous(labels = abs)+
#   labs(x = "Count", y = NULL, fill = "Comparison") +  
#   theme_minimal() + 
#   ggtitle("Number of changed TAD") +
#   scale_fill_manual(values = c("Down" = "skyblue", "Up" =  "salmon")) +  
#   theme(  
#     axis.title.x = element_text(size = 14),     
#     axis.title.y = element_text(size = 14),    
#     axis.text.x = element_text(hjust = 1, size = 12),   
#     axis.text.y = element_text(size = 12),    
#     plot.title = element_text(size = 16, face = "bold"),
#     legend.position = "bottom"
#   ) 
