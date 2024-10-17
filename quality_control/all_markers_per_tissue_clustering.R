rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(patchwork)
library(edgeR)
library(MASS) 
library(Seurat)
library(gridExtra)
library(ggrepel)
library(stringr)
tissue <- "iWAT"
CUTTag_search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
ATAC_search_table <- read.csv("data/samples/all/ATAC_search_table.csv")
search_table <- rbind(CUTTag_search_table,ATAC_search_table)
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}

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
    }
  }
  return(tissue_label)
}
all_markers_pca <- function(tissue){
  tab = read.delim(paste0("data/samples/",tissue,"/all_antibodys_10kb_bins.counts"),skip=1)  
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames <- colnames(tab)[7:length(tab)]
  new_colnames <- gsub(pattern, "\\1", colnames)
  colnames(tab)[7:length(tab)] <- new_colnames
  sorted_index <- order(new_colnames)
  order_colnames <- new_colnames[sorted_index] 
  counts <- tab[,order_colnames] 
  y= DGEList(counts=counts)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  logCPMs <- cpm(y, log = TRUE)
  pca <- prcomp(t(logCPMs))
  to_plot <- data.frame(pca$x)
  to_plot$sample_name <- rownames(to_plot)
  to_plot <- merge(to_plot, search_table, by="sample_name")
  to_plot$label <- paste0(to_plot$sample_name,"-",to_plot$mouse_ID,"-",to_plot$age)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  use.pcs <- c(1,2)
  labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))
  color <- read.table("data/samples/7_distinct_color.txt")
  color <- setNames(color$V1,c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3","ATAC"))
  to_plot$age <- factor(to_plot$age, levels = c("3m","24m"))
  p <- ggplot(to_plot, aes(x=PC1, y=PC2, color=antibody,shape=age)) + 
    geom_point(size=5) +theme_bw()+ scale_color_manual(values = color)+
    xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
    geom_text_repel(  
      data = to_plot,  
      aes(x = PC1, y = PC2, label = label, color = antibody),  
      size = 3,  
      box.padding = unit(0.35, "lines"),  
      point.padding = unit(0.3, "lines")  
    ) +
    ggtitle(paste(tissue_label_change(tissue)))
  return(p)
}

all_markers_correlation <- function(tissue){
  tab = read.delim(paste0("data/samples/",tissue,"/all_antibodys_10kb_bins.counts"),skip=1)  
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames <- colnames(tab)[7:length(tab)]
  new_colnames <- gsub(pattern, "\\1", colnames)
  colnames(tab)[7:length(tab)] <- new_colnames
  sorted_index <- order(new_colnames)
  order_colnames <- new_colnames[sorted_index] 
  counts <- tab[,order_colnames] 
  y= DGEList(counts=counts)
  keep = which(rowSums(cpm(y)>1)>=2)
  y = y[keep,]
  logCPMs <- cpm(y, log = TRUE)
  correlation_matrix <- cor(logCPMs,method = "spearman")
  annotation <- search_table[which(search_table$sample_name %in% colnames(correlation_matrix)),]
  rownames(annotation) <- annotation$sample_name
  annotation <- annotation[,c("antibody","age")]
  annotation$age <- factor(annotation$age,c("3m","24m"))
  color <- read.table("data/samples/7_distinct_color.txt")
  antibody_color <- setNames(color$V1,c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3","ATAC"))
  age_color <- setNames(c("red","blue"),c("3m","24m"))
  annotation_colors <- list(age=age_color,antibody=antibody_color)
  pheatmap::pheatmap(correlation_matrix,annotation = annotation, 
                     annotation_colors = annotation_colors, main = tissue_label_change(tissue))
  
}

all_markers_umap <- function(tissue){
  tab = read.delim(paste0("data/samples/",tissue,"/all_antibodys_10kb_bins.counts"),skip=1)  
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
  colnames <- colnames(tab)[7:length(tab)]
  new_colnames <- gsub(pattern, "\\1", colnames)
  colnames(tab)[7:length(tab)] <- new_colnames
  sorted_index <- order(new_colnames)
  order_colnames <- new_colnames[sorted_index] 
  counts <- tab[,order_colnames] 
  annotation <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  rownames(annotation) <- annotation$sample_name
  annotation <- annotation[,c("antibody","age","mouse_ID")]
  pbmc <- CreateSeuratObject(count=counts)
  pbmc <- AddMetaData(pbmc, metadata = annotation)  
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 10000)
  pbmc <- ScaleData(pbmc)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 10)
  pbmc <- RunUMAP(pbmc, dims = 1:10,n.neighbors = 4)
  pbmc <- FindNeighbors(pbmc, dims = 1:10)
  pbmc <- FindClusters(pbmc, resolution = 0.5,verbose = FALSE)
  pbmc$label <- paste0(rownames(pbmc@meta.data),"-", pbmc$mouse_ID,"-",pbmc$age)
  color <- read.table("data/samples/7_distinct_color.txt")
  antibody_color <- setNames(color$V1,c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3","ATAC"))
  cells <- as.data.frame(Embeddings(pbmc, "umap"))  
  cells$antibody <- pbmc@meta.data$antibody  
  cells$age <- pbmc@meta.data$age  
  cells$mouse_ID <- pbmc@meta.data$mouse_ID
  cells$label<- paste0(rownames(cells),"-", cells$mouse_ID,"-",cells$age)
  pbmc$age <- factor(pbmc$age,levels=c("3m","24m"))
  p <- DimPlot(pbmc,group.by = "antibody",label = F, shape.by = "age",pt.size = 3) +
    scale_color_manual(values = antibody_color) + 
    ggtitle(tissue_label_change(tissue))+
    geom_text_repel(data = cells, aes(x = umap_1, y = umap_2, label = label), size = 2.5)  
  return(p)
}

tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder","CB","jejunum","uterus","ovary")
p_list <- list()
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  p_list[[i]] <- all_markers_umap(tissue)
}
combined_plot <- plot_a_list(p_list,4,6)
ggsave("result/all/pca/all_tissues_plot/all_markers_per_tissue_umap.png",combined_plot,width = 40,height = 20,type="cairo")

tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder","CB","jejunum","uterus","ovary")
p_list <- list()
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  p_list[[i]] <- all_markers_pca(tissue)
}
combined_plot <- plot_a_list(p_list,4,6)
ggsave("result/all/pca/all_tissues_plot/all_markers_per_tissue_PCA.png",combined_plot,width = 40,height = 20,type="cairo")
