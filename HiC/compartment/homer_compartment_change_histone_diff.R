rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
library(grid)
tissue <- "lung"
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
diff_analysis <- function(tissue,antibody){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_compartment_50000.counts"),skip=1)
  counts = tab[,c(7:ncol(tab))]
  rownames(counts)= tab$Geneid
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  counts <- counts[,search_table$sample_name]
  search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
  search_table$age <- factor(search_table$age, levels = c("3m","24m"))
  
  age <- as.character(search_table$age)
  batch <- as.character(search_table$batch)
  mouse_ID <- search_table$mouse_ID
  age[which(age=="3m")] <- "young"
  age[which(age=="24m")] <- "old"
  colnames(counts) <- paste0(colnames(counts),"-",age,"-",mouse_ID,"-",batch)
  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  y$samples$year <- age
  y$samples$year <- factor(y$samples$year,c("young","old"))
  y$samples$batch <- search_table$batch
  
  y <- calcNormFactors(y)
  design <- model.matrix(~batch+year, y$samples)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = which(colnames(design) == "yearold"))
  tab<-tab[keep,]
  
  out = cbind(tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
              "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
              "LogFC.old-young"=lrt$table$logFC)
  
  out$Significant <- ifelse(out$`FDR.old-young` < 0.05 & abs(out$`LogFC.old-young`) >= log2(1.2), 
                            ifelse(out$`LogFC.old-young` > log2(1.2), "Up", "Down"), "Stable")
  
  write.csv(out,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_diff_in_compartment_50000.csv"))
  
  compartment <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_50000.csv"))
  compartment$Geneid <- paste0(compartment$chr,":",compartment$start,"-",compartment$end)
  compartment <- compartment[,c("Geneid","condition")]
  out <- out[,c("LogFC.old-young","logCPM","Significant","Geneid")]
  
  df <- merge(compartment,out,by="Geneid")
  to_plot <- df 
  to_plot$condition <- factor(to_plot$condition, levels=c("A-A","B-B","A-B","B-A"))
  p <-  ggplot(to_plot, aes(x = condition , y = `LogFC.old-young`, fill= condition)) +  
    geom_boxplot() +
    theme_minimal()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ylab(paste0(antibody," log2(Fold change)")) +
    xlab("H3K9me3 changes")+
    ggtitle(tissue_label_change(tissue))
  ggsave(paste0("result/all/diff/",antibody,"/",tissue,"_relationship_change_in_compartment_with_compartment_change.png"),p,width=4,height=5,type="cairo")
}
tissues <- c("brain","CB", "kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip", "mammarygland", "stomach", "thymus","skin","muscle")
antibodys <- c("H3K9me3","H3K36me3","H3K27me3","H3K4me3","H3K4me1","H3K27ac")
for(tissue in tissues){
  for(antibody in antibodys){
    diff_analysis(tissue,antibody)
  }
}
for(antibody in antibodys){
  summary <- data.frame() 
  for(tissue in tissues){
    active_mark <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_diff_in_compartment_50000.csv"),row.names = 1)
    active_mark <- active_mark[,c("Geneid","LogFC.old.young","logCPM","Significant")]
    colnames(active_mark)[4] <- "active_mark_significant"
    compartment <- read.csv(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/compartment_change_50000.csv"))
    compartment$Geneid <- paste0(compartment$chr,":",compartment$start,"-",compartment$end)
    compartment <- compartment[,c("Geneid","condition")]
    df <- merge(compartment,active_mark,by="Geneid")
    to_plot <- df 
    to_plot$condition <- factor(to_plot$condition, levels=c("A-A","A-B","B-B","B-A"))
    to_plot$tissue <- tissue_label_change(tissue)
    summary <- rbind(summary,to_plot)
  }
  to_plot <- summary %>%  
    group_by(tissue, condition) %>%  
    summarise(  
      median_logFC = median(LogFC.old.young, na.rm = TRUE),  
      count = n()  
    ) %>%   
    mutate(median_logFC = ifelse(count < 10, NA, median_logFC))
  
  p_value_summary <- data.frame(
    'A-A' = rep(NA, length(tissues)), 
    'A-B' = rep(NA, length(tissues)), 
    'B-B' = rep(NA, length(tissues)),
    'B-A' = rep(NA, length(tissues))
  )
  colnames(p_value_summary) <- c("A-A","A-B","B-B","B-A")
  rownames(p_value_summary) <- sapply(tissues, tissue_label_change)
  for(tissue in tissues){
    t_A_A_summary <-  summary[which(summary$tissue == tissue_label_change(tissue) & summary$condition=="A-A"),]
    t_B_B_summary <- summary[which(summary$tissue == tissue_label_change(tissue) & summary$condition=="B-B"),]
    t_A_B_summary <- summary[which(summary$tissue == tissue_label_change(tissue) & summary$condition=="A-B"),]
    t_B_A_summary <- summary[which(summary$tissue == tissue_label_change(tissue) & summary$condition=="B-A"),]
    
    if(nrow(t_A_B_summary) > 10){
      test <- wilcox.test(t_A_B_summary$LogFC.old.young, t_A_A_summary$LogFC.old.young)
      p_value_summary[tissue_label_change(tissue),"A-B"] <- test$p.value
    }
    if(nrow(t_B_A_summary) > 10){
      test <- wilcox.test(t_B_A_summary$LogFC.old.young,t_B_B_summary$LogFC.old.young)
      p_value_summary[tissue_label_change(tissue),"B-A"] <- test$p.value
    }
  }
  mark_significance <- function(p_value) {
    if (is.na(p_value)) {
      return(NA)
    } else if (p_value < 0.001) {
      return("***")
    } else if (p_value < 0.01) {
      return("**")
    } else if (p_value < 0.05) {
      return("*")
    } else {
      return(NA)
    }
  }
  
  p_value_summary <- p_value_summary %>%
    mutate(
      `A-A` = sapply(`A-A`, mark_significance),
      `A-B` = sapply(`A-B`, mark_significance),
      `B-B` = sapply(`B-B`, mark_significance),
      `B-A` = sapply(`B-A`, mark_significance)
    )
  p_value_summary$tissue <- rownames(p_value_summary)
  p_value_long <- reshape2::melt(p_value_summary,id.vars = "tissue")
  names(p_value_long) <- c("Tissue", "Type", "Label")
  colnames(to_plot)[c(1:3)] <- c("Tissue","Type","Value")
  merged_data <- merge(to_plot, p_value_long, by = c("Tissue", "Type"), all.x = TRUE)
  merged_data$Value[which(merged_data$Value > 1)] <- 1
  merged_data$Value[which(merged_data$Value < -1)] <- -1
  
  ggplot(merged_data, aes(x = Type, y = Tissue, fill = Value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         limits = c(-1, 1), midpoint = 0) +
    theme_minimal() +
    ggtitle(antibody)+
    geom_text(aes(label = Label), color = "black", size = 4, na.rm = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
