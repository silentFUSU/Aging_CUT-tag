rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
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
      tissue_label <- "IWAT"
    }
  }
  return(tissue_label)
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
get_all_tissues_bcv <- function(tissues,antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  bcv_summary <- data.frame()
  for(tissue in tissues){
    search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
    tab=read.delim(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),skip=1)
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
    keep = which(rowSums(edgeR::cpm(y)>1)>=2)
    y = y[keep,]
    y$samples$year <- age
    y$samples$year <- factor(y$samples$year,c("young","old"))
    y$samples$batch <- search_table$batch
    
    y <- calcNormFactors(y)
    design <- model.matrix(~year+batch, y$samples)
    y<-estimateCommonDisp(y)
    y<-estimateGLMTagwiseDisp(y,design)
    bcv <- data.frame(peaks=rownames(counts[keep,]),
                      bcv=sqrt(y$tagwise.dispersion))
    bcv$tissue <- tissue_label_change(tissue)
    bcv_summary <- rbind(bcv_summary,bcv)
  }
  color <- read.table("data/samples/30_distinct_color.txt")
  color <- color$V1
  color <- setNames(color,sort(unique(bcv_summary$tissue)))
  p <- ggplot(bcv_summary,aes(x=tissue,y=bcv,fill=tissue))+
    geom_violin()+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    scale_fill_manual(values = color) +
    theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+ylab("bcv")+
    xlab("")+labs(fill = "", color = "")+ggtitle(paste0("bcv each tissue ",antibody))+ylim(0,1)
  ggsave(paste0("result/all/QC/bcv/all_tissues_",antibody,"_bcv_remove_batch_effect.png"),p,width = 18,height = 10,type="cairo")  
  bcv_average <- mean(bcv_summary$bcv)
  saveRDS(bcv_summary,"data/samples/all/bcv.rds")
  }

color <- read.table("data/samples/30_distinct_color.txt")
color <- color$V1
# bcv <- readRDS("data/samples/all/bcv.rds")
bcv <- bcv_summary
bcv$tissue[which(bcv$tissue=="IWAT")] <- "iWAT"
color <- setNames(color,sort(unique(bcv$tissue)))
# bcv$tissue <- factor(bcv$tissue,levels = c("Testis","Tongue","Stomach","Cecum","Colon","Pancreas","Ileum",
#                                            "Liver","Heart","Jejunum","Hippocampus","Muscle","Bone Marrow",
#                                            "Ovary","iWAT","Cortex","Aorta","Uterus","Bladder","Spleen",
#                                            "Thymus","Kidney","Skin","Cerebellum","Lung","BAT","Mammary Gland"))
bcv$tissue <- factor(bcv$tissue,levels = c("Pancreas","Cecum","Colon","Kidney","Ileum","Spleen","Testis","Stomach",
                                           "Skin","Jejunum","Tongue","iWAT","Liver","Bone Marrow","Heart","Cerebellum",
                                           "Uterus","Lung","Thymus","Bladder","BAT","Hippocampus","Aorta","Muscle",
                                           "Ovary","Cortex","Mammary Gland"))
ggplot(bcv,aes(x=bcv,y=tissue,fill=tissue))+
  geom_violin()+
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
  scale_fill_manual(values = color) +
  theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust = 1))+ylab("")+
  xlab("bcv ~age+batch")+labs(fill = "", color = "")+ggtitle(paste0("bcv each tissue ",antibody))+xlim(0,1)




