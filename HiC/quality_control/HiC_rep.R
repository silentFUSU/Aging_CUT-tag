rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(hicrep)
library(stringr)
library(strawr)
options(bitmapType = "cairo")
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
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
# tissue <- "lung"
search_table <- read.csv("data/samples/all/HiC_search_table.csv")
# search_table <- search_table[which(search_table$tissue==tissue),]
cor_df <- as.data.frame(matrix(data = NA, nrow = nrow(search_table), ncol = nrow(search_table)))
colnames(cor_df) <- search_table$sample_name
rownames(cor_df) <- search_table$sample_name
for(i in c(1:nrow(search_table))){
  for(j in c(1:nrow(search_table))){
    scc <- 0
    for (chr in paste0("chr", c(as.character(1:19), "X","Y"))){
      mat1 <- hic2mat(paste0("data/samples/HiC/",search_table[i,"tissue"],"/juicer/",search_table[i,"sample_name"],".allValidPairs.hic"), chromosome1 = chr, chromosome2 = chr, resol =   100000, method = "NONE") 
      mat2 <- hic2mat(paste0("data/samples/HiC/",search_table[j,"tissue"],"/juicer/",search_table[j,"sample_name"],".allValidPairs.hic"), chromosome1 = chr, chromosome2 = chr, resol =   100000, method = "NONE") 
      scc.out = get.scc(mat1, mat2, resol = 100000, h = 5, lbr = 0, ubr = 5000000)
      scc <- scc + scc.out$scc
    }
    scc <- scc/length(c(as.character(1:19), "X","Y"))
    cor_df[search_table[i,"sample_name"],search_table[j,"sample_name"]] <- scc
  }
}
annotation <- search_table[,c("sample_name","age","tissue")]
annotation$age[which(annotation$age=="3M")] <- "Young"
annotation$age[which(annotation$age=="24M")] <- "Old"
annotation$age <- factor(annotation$age,levels=c("Young","Old"))
rownames(annotation) <- annotation$sample_name
annotation <- annotation[,-1]
annotation$tissue <- sapply(annotation$tissue, tissue_label_change) 
colors <- read.table("data/samples/20_distinct_color.txt")
annotation_color <- list(tissue=setNames(colors$V1[1:length(unique(annotation$tissue))],unique(annotation$tissue)))
pheatmap::pheatmap(cor_df,annotation_row =annotation,annotation_colors = annotation_color,main = paste0("HiC replicate scc"),filename = paste0("result/HiC/all_tissues_HiCrep_scc_heatmap.png"),width = 12,height = 10)


tissue <- "Hip"
per_tissue_scc <- function(tissue){
  search_table <- read.csv("data/samples/all/HiC_search_table.csv")
  search_table <- search_table[which(search_table$tissue==tissue),]
  cor_df <- as.data.frame(matrix(data = NA, nrow = nrow(search_table), ncol = nrow(search_table)))
  colnames(cor_df) <- search_table$sample_name
  rownames(cor_df) <- search_table$sample_name
  for(i in c(1:nrow(search_table))){
    for(j in c(1:nrow(search_table))){
      scc <- 0
      for (chr in paste0("chr", c(as.character(1:19), "X","Y"))){
        mat1 <- hic2mat(paste0("data/samples/HiC/",search_table[i,"tissue"],"/juicer/",search_table[i,"sample_name"],".allValidPairs.hic"), chromosome1 = chr, chromosome2 = chr, resol =   100000, method = "NONE") 
        mat2 <- hic2mat(paste0("data/samples/HiC/",search_table[j,"tissue"],"/juicer/",search_table[j,"sample_name"],".allValidPairs.hic"), chromosome1 = chr, chromosome2 = chr, resol =   100000, method = "NONE") 
        scc.out = get.scc(mat1, mat2, resol = 100000, h = 5, lbr = 0, ubr = 5000000)
        scc <- scc + scc.out$scc
      }
      scc <- scc/length(c(as.character(1:19), "X","Y"))
      cor_df[search_table[i,"sample_name"],search_table[j,"sample_name"]] <- scc
    }
  }
  annotation <- search_table[,c("sample_name","age")]
  rownames(annotation) <- annotation$sample_name
  annotation$age[which(annotation$age=="3M")] <- "Young"
  annotation$age[which(annotation$age=="24M")] <- "Old"
  annotation$age <- factor(annotation$age,levels=c("Young","Old"))
  annotation <- annotation[,-1,drop=F]
  dir.create(paste0("result/HiC/",tissue))
  dir.create(paste0("result/HiC/",tissue,"/quality_control/"))
  pheatmap::pheatmap(cor_df,annotation_row = annotation,main = paste0(tissue_label_change(tissue)),filename = paste0("result/HiC/",tissue,"/quality_control/HiCrep_scc_heatmap.png"),width = 7,height = 6)
}
tissues <- c("lung","liver","thymus","CB","brain","colon","stomach","heart","kidney","bonemarrow")
for(tissue in tissues){
  per_tissue_scc(tissue)  
}
