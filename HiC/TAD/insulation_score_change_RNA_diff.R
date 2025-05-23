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
  rownames(tab) <- tab$Geneid
  tab <- tab[,-1]
  colnames <- colnames(tab)[6:length(tab)]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+).*"
  colnames(tab)[6:length(tab)] <- gsub(pattern, "\\1", colnames(tab)[6:length(tab)] )
  counts <- tab[6:length(tab)]
  group <- read.csv("data/samples/RNA/sample_tissue_info.csv",sep = ',')
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  colnames(group)[1] <- "sample_name"
  group <- merge(group,search_table,by="sample_name")
  group <- group[which(group$sample_name %in% colnames(counts)),]
  group$sample_name <- factor(group$sample_name,levels = colnames(counts))
  group <- group[order(group$sample_name),]
  age <- group$Age
  age[which(age=="3m")] <- "young"
  age[which(age=="24m")] <- "old"
  mouse_ID <- group$mouse_ID
  
  colnames(counts) <- paste0(colnames(counts),"-",mouse_ID,"-",age)
  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
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
  
  tad <- read.csv(paste0("data/samples/HiC/",tissue,"/TAD/insulation_score/",tissue,"_redundant_20000_TAD_diff.csv"))
  tad <- tad %>%
    separate(X, into = c("chr", "start", "end"), sep = "-", convert = TRUE)
  tad$start <- tad$start-1
  tad$Geneid <- paste0(tad$chr,":",tad$start,"-",tad$end)
  tad <- tad[,c("Geneid","Significant")]
  out$Geneid <- rownames(out)
  out <- out[,c("logFC","logCPM","Significant","Geneid")]
  
  df <- merge(tad,out,by="Geneid")
  to_plot <- df 
  to_plot$condition <- factor(to_plot$condition, levels=c("A-A","B-B","A-B","B-A"))
  p <- ggplot(to_plot, aes(x = condition , y = logFC, fill=condition)) +  
    geom_boxplot() +
    theme_minimal()+
    theme(text = element_text(size = 20),legend.position = "none")+
    ylab("RNA log2(Fold change)") +
    xlab("compartment changes")+
    ggtitle(tissue_label_change(tissue))
  ggsave(paste0("result/RNA/",tissue,"/relationship_gene_expression_in_compartment_with_compartment_change.png"),p,width=4,height=5,type="cairo")
}