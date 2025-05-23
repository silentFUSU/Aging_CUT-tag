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
diff_expression_analysis <- function(tissue,antibody){
  search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
  tab = read.delim(paste0("data/samples/",tissue,"/",antibody,"/",tissue,"_H3K9me3_peaks.counts"),skip=1)
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
  
  write.csv(out,paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_diff_in_H3K9me3_peaks_remove_batch_effect.csv"))
  
  histone <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv"))
  histone <- histone[,c("Geneid","LogFC.old.young","Significant")]
  out <- out[,c("LogFC.old-young","logCPM","Significant","Geneid")]
  
  df <- merge(histone,out,by="Geneid")
  to_plot <- df 
  
  colnames(to_plot) <- c("Geneid","H3K9me3_logFC","Significant","active_logFC","active_logCPM","active_Significant")
  # to_plot <- to_plot[which(to_plot$Significant != "Stable"),]
  to_plot$Significant <- factor(to_plot$Significant, levels=c("Up","Stable","Down"))
  p <- ggplot(to_plot, aes(x =Significant , y = active_logFC, fill=Significant)) +  
    geom_boxplot() +
    theme_minimal()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ylab(paste0(antibody," log2(Fold change)")) +
    xlab("H3K9me3 changes")+
    ggtitle(tissue_label_change(tissue))
  ggsave(paste0("result/all/diff/",antibody,"/",tissue,"_relationship_change_in_H3K9me3_peaks_with_H3K9me3_changes.png"),p,width=4,height=5,type="cairo")
}
tissues <- sort(c("BAT","mammarygland","CB","lung","kidney","aorta","brain","spleen",
                  "thymus","skin","bladder","bonemarrow","Hip","heart",
                  "muscle","jejunum","uterus","ovary","liver","tongue",
                  "cecum","colon","testis","stomach","pancreas","iWAT","ileum"))
antibody <- "H3K36me3"
for(tissue in tissues){
  diff_expression_analysis(tissue,antibody)
}

summary <- data.frame()
for(tissue in tissues){
  active_mark <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_diff_in_H3K9me3_peaks_remove_batch_effect.csv"),row.names = 1)
  active_mark <- active_mark[,c("Geneid","LogFC.old.young","logCPM","Significant")]
  colnames(active_mark)[4] <- "active_mark_significant"
  histone <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W5000-G10000-E100_recursion_diff_after_remove_batch_effect.csv"))
  histone <- histone[,c("Geneid","Length","LogFC.old.young","Significant")]
  colnames(histone)[4] <- "histone_significant"
  to_plot <- merge(histone,active_mark,by="Geneid")
  colnames(to_plot)[c(3,5)] <- c("LogFC.old.young","logFC")
  to_plot$tissue <- tissue_label_change(tissue)
  summary <- rbind(summary,to_plot)
}

color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
color <- setNames(color,sort(unique(summary$tissue)))
summary <- summary[which(summary$Length > 200000),]
summary$tissue[which(summary$tissue=="IWAT")] <- "iWAT"
tissues_order <- c("Testis","Tongue","Stomach","Cecum","Colon","Pancreas","Ileum","Liver","Heart","Jejunum","Hippocampus","Muscle","Bone Marrow","Ovary","iWAT","Cortex","Aorta","Uterus","Bladder","Spleen","Thymus","Kidney","Skin","Cerebellum","Lung","BAT","Mammary Gland")

decrease_summary <- summary[which(summary$histone_significant=="Down"),]
filter <- as.data.frame(table(decrease_summary$tissue))
filter <- filter$Var1[which(filter$Freq > 10)]
decrease_summary$tissue <- factor(decrease_summary$tissue,levels=tissues_order)
decrease_summary$logFC[which(!decrease_summary$tissue %in% filter)] <- NA
if(length(tissues_order[-which(tissues_order %in% decrease_summary$tissue)])>0){
  for(tissue in tissues_order[-which(tissues_order %in% decrease_summary$tissue)]){
    decrease_summary <- rbind(decrease_summary, data.frame(Geneid=NA,
                                                           Length=NA,
                                                           LogFC.old.young=NA,
                                                           histone_significant=NA,
                                                           logFC=NA,
                                                           logCPM=NA,
                                                           active_mark_significant=NA,
                                                           tissue=tissue))
  }
}

ggplot(decrease_summary, aes(x = tissue, y = logFC, fill=tissue)) +  
  geom_boxplot(outlier.shape = NA) +
  theme_minimal()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme(  
    text = element_text(size = 20),  
    legend.position = "none",  
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  scale_fill_manual(values = color) +
  ylab(paste0(antibody," log2(Fold change)")) +
  xlab(NULL)+
  ggtitle("Regions with decreased H3K9me3")+
  ylim(-3,3)

Stable_summary <- summary[which(summary$histone_significant=="Stable"),]
filter <- as.data.frame(table(Stable_summary$tissue))
filter <- filter$Var1[which(filter$Freq > 10)]
Stable_summary$tissue <- factor(Stable_summary$tissue,levels=tissues_order)
Stable_summary$logFC[which(!Stable_summary$tissue %in% filter)] <- NA
if(length(tissues_order[-which(tissues_order %in% Stable_summary$tissue)])>0){
  for(tissue in tissues_order[-which(tissues_order %in% Stable_summary$tissue)]){
    Stable_summary <- rbind(Stable_summary, data.frame(Geneid=NA,
                                                       Length=NA,
                                                       LogFC.old.young=NA,
                                                       histone_significant=NA,
                                                       logFC=NA,
                                                       logCPM=NA,
                                                       active_mark_significant=NA,
                                                       tissue=tissue))
  }
}

ggplot(Stable_summary, aes(x = tissue, y = logFC, fill=tissue)) +  
  geom_boxplot(outlier.shape = NA) +
  theme_minimal()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme(  
    text = element_text(size = 20),  
    legend.position = "none",  
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  scale_fill_manual(values = color) +
  ylab(paste0(antibody," log2(Fold change)")) +
  xlab(NULL)+
  ggtitle("Regions with Stable H3K9me3")+
  ylim(-3,3)

increase_summary <- summary[which(summary$histone_significant=="Up"),]
filter <- as.data.frame(table(increase_summary$tissue))
filter <- filter$Var1[which(filter$Freq > 10)]
increase_summary$tissue <- factor(increase_summary$tissue,levels=tissues_order)
increase_summary$logFC[which(!increase_summary$tissue %in% filter)] <- NA
if(length(tissues_order[-which(tissues_order %in% increase_summary$tissue)])>0){
  for(tissue in tissues_order[-which(tissues_order %in% increase_summary$tissue)]){
    increase_summary <- rbind(increase_summary, data.frame(Geneid=NA,
                                                           Length=NA,
                                                           LogFC.old.young=NA,
                                                           histone_significant=NA,
                                                           logFC=NA,
                                                           logCPM=NA,
                                                           active_mark_significant=NA,
                                                           tissue=tissue))
  }
}
ggplot(increase_summary, aes(x = tissue, y = logFC, fill=tissue)) +  
  geom_boxplot(outlier.shape = NA) +
  theme_minimal()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme(  
    text = element_text(size = 20),  
    legend.position = "none",  
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  scale_fill_manual(values = color) +
  ylab(paste0(antibody," log2(Fold change)")) +
  xlab(NULL)+
  ylim(-3,3)+
  ggtitle("Regions with Increase H3K9me3")

to_plot <- summary %>%  
  group_by(tissue, histone_significant) %>%  
  summarise(  
    median_logFC = median(logFC, na.rm = TRUE),  
    count = n()  
  ) %>%   
  mutate(median_logFC = ifelse(count < 10, NA, median_logFC))
to_plot <- to_plot[,-ncol(to_plot)]
to_plot <- as.data.frame(to_plot)
to_plot <- reshape2::dcast(to_plot, tissue ~ histone_significant, value.var = "median_logFC")  
rownames(to_plot) <- to_plot$tissue
to_plot <- to_plot[,-1]
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  
breaks <- c(seq(-0.5, -0.21, length.out = 40), seq(-0.2, 0.2, length.out = 20), seq(0.21, 0.5, length.out = 40))  
to_plot <- to_plot[tissues_order,]
to_plot <- to_plot[,c("Up","Stable","Down")]
pheatmap::pheatmap(to_plot,breaks = breaks,color = color_palette,cluster_rows = F,cluster_cols = F,na_col = "grey")

p_value_summary <- data.frame(
  Up = rep(NA, 27), 
  Stable = rep(NA, 27), 
  Down = rep(NA, 27)  
)
rownames(p_value_summary) <- sapply(tissues, tissue_label_change)
for(tissue in tissues){
  t_Up_summary <-  summary[which(summary$tissue == tissue_label_change(tissue) & summary$histone_significant=="Up"),]
  t_Down_summary <- summary[which(summary$tissue == tissue_label_change(tissue) & summary$histone_significant=="Down"),]
  t_Stable_summary <- summary[which(summary$tissue == tissue_label_change(tissue) & summary$histone_significant=="Stable"),]
  if(nrow(t_Up_summary) > 10){
    test <- wilcox.test(t_Up_summary$logFC,t_Stable_summary$logFC)
    p_value_summary[tissue_label_change(tissue),"Up"] <- test$p.value
  }
  if(nrow(t_Down_summary) > 10){
    test <- wilcox.test(t_Down_summary$logFC,t_Stable_summary$logFC)
    p_value_summary[tissue_label_change(tissue),"Down"] <- test$p.value
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
    Up = sapply(Up, mark_significance),
    Stable = sapply(Stable, mark_significance),
    Down = sapply(Down, mark_significance)
  )
p_value_summary <- p_value_summary[tissues_order,]
to_plot$tissue <- rownames(to_plot)
df_long <- reshape2::melt(to_plot)
names(df_long) <- c("Tissue", "Type", "Value")

p_value_summary$tissue <- rownames(p_value_summary)
p_value_long <- reshape2::melt(p_value_summary,id.vars = "tissue")
names(p_value_long) <- c("Tissue", "Type", "Label")
merged_data <- merge(df_long, p_value_long, by = c("Tissue", "Type"), all.x = TRUE)
merged_data$Tissue <- factor(merged_data$Tissue,levels=rev(tissues_order))
ggplot(merged_data, aes(x = Type, y = Tissue, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme_minimal() +
  xlab("H3K9me3 change condition")+
  ggtitle(paste0(antibody," log2(Fold change)"))+
  geom_text(aes(label = Label), color = "black", size = 4, na.rm = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
