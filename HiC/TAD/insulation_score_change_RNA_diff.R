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

diff_expression_analysis <- function(tissue){
  tab = read.delim(paste0("data/samples/RNA/",tissue,"/counts/",tissue,"_20000_redundant_tads.counts"),skip=1)
  if(tissue %in% c("mammarygland","ovary","uterus")){
    tab <- tab[which(tab$Chr %in% paste0("chr",c(1:19,"X"))),]
  }
  rownames(tab) <- tab$Geneid
  tab <- tab[,-1]
  colnames <- colnames(tab)[6:length(tab)]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|HM[0-9]+).*"
  colnames(tab)[6:length(tab)] <- gsub(pattern, "\\1", colnames(tab)[6:length(tab)] )
  counts <- tab[6:length(tab)]
  group <- read.csv("data/samples/RNA/sample_tissue_info.csv",sep = ',')
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  colnames(group)[1] <- "sample_name"
  group <- merge(group,search_table,by="sample_name")
  group <- group[which(group$sample_name %in% colnames(counts)),]
  counts <- counts[,which(colnames(counts) %in% group$sample_name)]
  group$sample_name <- factor(group$sample_name,levels = colnames(counts))
  group <- group[order(group$sample_name),]
  age <- group$Age
  age[which(age=="3m")] <- "young"
  age[which(age=="24m")] <- "old"
  mouse_ID <- group$mouse_ID
  
  colnames(counts) <- paste0(colnames(counts),"-",mouse_ID,"-",age)
  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>0)>=2)
  y = y[keep,]
  y$samples$group <- factor(y$samples$group, levels=c("young","old"))
  design <- model.matrix(~group, y$samples)
  y <- calcNormFactors(y)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = 2)
  tab<-tab[keep,]
  out = cbind(cpm(y),lrt$table, "fdr"=p.adjust(lrt$table$PValue,method="BH"))
  out$Significant <- ifelse(out$fdr< 0.05 & abs(out$logFC) >= 0, 
                            ifelse(out$logFC > 0, "Up", "Down"), "Stable")
  write.csv(out,paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_in_20000_redundant_tads.csv"))
  
  tad <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_20000_TAD_diff_larger_250000.csv"))
  tad <- tad %>%
    separate(X, into = c("chr", "start", "end"), sep = "-", convert = TRUE)
  tad$start <- tad$start-1
  tad$Geneid <- paste0(tad$chr,":",tad$start,"-",tad$end)
  tad <- tad[,c("Geneid","Significant")]
  out$Geneid <- rownames(out)
  out <- out[,c("logFC","logCPM","Geneid")]
  
  df <- merge(tad,out,by="Geneid")
  to_plot <- df 
  to_plot$condition <- factor(to_plot$Significant, levels=c("Up","Stable","Down"))
  p <- ggplot(to_plot, aes(x = condition , y = logFC, fill=condition)) +  
    geom_boxplot() +
    theme_minimal()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ylab("RNA log2(Fold change)") +
    xlab("compartment changes")+
    ggtitle(tissue_label_change(tissue))
  ggsave(paste0("result/RNA/",tissue,"/relationship_gene_expression_in_20000_redundant_tads_with_20000_redundant_tads_change.png"),p,width=4,height=5,type="cairo")
}
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","skin","muscle","cecum","ileum","pancreas","spleen")
for(tissue in tissues){
  print(tissue)
  diff_expression_analysis(tissue)
}


summary <- data.frame()
for(tissue in tissues){
  gene <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_change_in_20000_redundant_tads.csv"))
  colnames(gene)[1] <- "Geneid"
  gene <- gene[,c("Geneid","logFC","logCPM","Significant")]
  colnames(gene)[4] <- "RNA_significant"
  tad <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_20000_TAD_diff_larger_250000.csv"))
  tad <- tad %>%
    separate(X, into = c("chr", "start", "end"), sep = "-", convert = TRUE)
  tad$start <- tad$start-1
  tad$Geneid <- paste0(tad$chr,":",tad$start,"-",tad$end)
  tad <- tad[,c("Geneid","Significant")]
  
  df <- merge(tad,gene,by="Geneid")
  to_plot <- df 
  to_plot$condition <- factor(to_plot$Significant, levels=c("Up","Stable","Down"))
  to_plot$tissue <- tissue_label_change(tissue)
  summary <- rbind(summary,to_plot)
}

to_plot <- summary %>%  
  group_by(tissue, condition) %>%  
  summarise(  
    median_logFC = median(logFC, na.rm = TRUE),  
    count = n()  
  ) %>%   
  mutate(median_logFC = ifelse(count < 20, NA, median_logFC))

p_value_summary <- data.frame(
  'Up' = rep(NA, length(tissues)), 
  'Stable' = rep(NA, length(tissues)), 
  'Down' = rep(NA, length(tissues))
)
rownames(p_value_summary) <- sapply(tissues, tissue_label_change)

for(tissue in tissues){
  t_Up_summary <-  summary[which(summary$tissue == tissue_label_change(tissue) & summary$condition=="Up"),]
  t_Stable_summary <- summary[which(summary$tissue == tissue_label_change(tissue) & summary$condition=="Stable"),]
  t_Down_summary <- summary[which(summary$tissue == tissue_label_change(tissue) & summary$condition=="Down"),]
  
  if(nrow(t_Up_summary) > 20){
    test <- wilcox.test(t_Up_summary$logFC, t_Stable_summary$logFC)
    p_value_summary[tissue_label_change(tissue),"Up"] <- test$p.value
  }
  if(nrow(t_Down_summary) > 20){
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
p_value_summary$tissue <- rownames(p_value_summary)
p_value_long <- reshape2::melt(p_value_summary,id.vars = "tissue")
names(p_value_long) <- c("Tissue", "Type", "Label")
colnames(to_plot)[c(1:3)] <- c("Tissue","Type","Value")
merged_data <- merge(to_plot, p_value_long, by = c("Tissue", "Type"), all.x = TRUE,all.y =TRUE)
ggplot(merged_data, aes(x = Type, y = Tissue, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       limits = c(-1, 1), midpoint = 0) +
  theme_minimal() +
  ggtitle("Gene expression")+
  geom_text(aes(label = Label), color = "black", size = 4, na.rm = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
