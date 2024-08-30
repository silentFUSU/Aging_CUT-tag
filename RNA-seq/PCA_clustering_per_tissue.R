rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
tissue<-"jejunum"
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
tissue_label_change <- function(tissue){
  if(tissue=="FC"){
    tissue_label <- "Cortex"
  }else if(tissue == "Hip"){
    tissue_label <- "Hippocampus"
  }else if(tissue == "CB"){
    tissue_label <- "Cerebellum"
  }else{
    tissue_label <- str_to_title(tissue)
    if(tissue_label == "Bonemarrow"){
      tissue_label <- "Bone Marrow"
    }
  }
  return(tissue_label)
} 

search_table <- read.csv("data/samples/all/RNA_search_table.csv")
PCA_per_tissue <- function(tissue){
  tab = read.delim(paste0("data/samples/RNA/",tissue,"/combined-chrM.nodup.counts"),skip=1)
  rownames(tab) <- tab$Geneid
  tab <- tab[,-1]
  colnames <- colnames(tab)[6:length(tab)]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+).*"
  new_colnames <- gsub(pattern, "\\1", colnames)
  colnames(tab)[6:length(tab)] <- new_colnames
  sorted_index <- order(new_colnames)
  order_colnames <- new_colnames[sorted_index] 
  counts <- tab[,order_colnames]   
  group <- read.csv("data/samples/RNA/sample_tissue_info.csv",sep = ',')
  group <- group[which(group$SampleID %in% colnames(counts)),]
  group <- group[order(group$SampleID),]
  age <- group[which(group$SampleID %in% colnames(counts)),"Age"]

  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  
  logCPMs <- cpm(y, log = TRUE)
  pca <- prcomp(t(logCPMs))
  to_plot <- data.frame(pca$x, age = paste0(y$samples$group))
  to_plot$rownames <- rownames(to_plot)
  table <- search_table[which(search_table$sample_name %in% to_plot$rownames),]
  to_plot$rownames <- paste0(table$sample_name,"-",table$mouse_ID,"-",table$age)
  to_plot$age <- factor(to_plot$age,levels=c("3m","24m"))
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  
 p <- ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
    geom_point(size=5) +theme_bw()+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
    geom_text_repel(  
      data = to_plot,  
      aes(x = PC1, y = PC2, label = rownames, color = age),  
      size = 5,  
      box.padding = unit(0.35, "lines"),  
      point.padding = unit(0.3, "lines")  
    ) +
    ggtitle(tissue_label_change(tissue))
 print(p)
  return(p)
}

tissues <- c("skin","CB","spleen","heart","bladder","tongue","uterus","aorta","thymus","stomach","Hip","FC","BAT","iWAT","muscle","bonemarrow","lung","kidney","liver","testis","colon","cecum","ileum","jejunum")
p_list <- list()
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  p_list[[i]] <- PCA_per_tissue(tissue)
}
combined_plot <- plot_a_list(p_list,4,6)
ggsave(paste0("result/RNA/per_tissue_PCA_nodup.png"),combined_plot,width = 35,height = 20,type="cairo")

