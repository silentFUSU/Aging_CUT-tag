rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(patchwork)
library(stringr)
library(dplyr)
library(tidyr)
library(UpSetR)
library(pheatmap)
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle")
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")){
    return ("10kb")
  }else{
    return("1kb")
  }  
}

jaccard_index<-data.frame(tissue = character(), 
                          jaccard_index = numeric(),
                          antibody1 = character(),
                          antibody1_condition=character(),
                          antibody2 = character(),
                          antibody2_condition=character(),
                          stringsAsFactors = FALSE)

conditions<-c("Up","Down")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")
for (i in c(1:length(tissues))){
  tissue<-tissues[i]
  for(j in c(1:length(antibodys))){
    antibody1 <- antibodys[j]
    for(k in c(1:length(antibodys))){
      antibody2 <- antibodys[k]
      df1 <- read.csv(paste0("data/samples/",tissue,"/",antibody1,"/",antibody1,"_",bin_size(antibody1),"_bins_diff_after_remove_batch_effect.csv")) 
      df2 <- read.csv(paste0("data/samples/",tissue,"/",antibody2,"/",antibody2,"_",bin_size(antibody2),"_bins_diff_after_remove_batch_effect.csv")) 
      for (condition1 in conditions){
        peaks1 <- df1$Geneid[which(df1$Significant_bar==condition1)]
        for (condition2 in conditions){
          peaks2 <- df2$Geneid[which(df2$Significant_bar==condition2)]
          t_jaccard_index <- data.frame(tissue = tissue, 
                                          jaccard_index = jaccard(peaks1,peaks2),
                                          antibody1 = antibody1,
                                          antibody1_condition=condition1,
                                          antibody2 = antibody2,
                                          antibody2_condition=condition2,
                                          stringsAsFactors = FALSE)
          jaccard_index <- rbind(jaccard_index,t_jaccard_index)
          }
        }
    }
  }
}
write.csv(jaccard_index,"result/all/jaccard_index.csv",row.names = F)
jaccard_index <-read.csv("result/all/jaccard_index.csv")
combinations <- paste(rep(antibodys, each = length(conditions)), rep(conditions, length(antibodys)), sep = "_")
mat <- matrix(0, nrow = length(combinations), ncol = length(combinations), dimnames = list(combinations, combinations))  
jaccard_index$label2 <- paste0(jaccard_index$antibody2,"_",jaccard_index$antibody2_condition)
jaccard_index$label1 <- paste0(jaccard_index$antibody1,"_",jaccard_index$antibody1_condition)
jaccard_index[is.na(jaccard_index)] <- 0


for(i in c(1:nrow(mat))){
  for(j in c(1:ncol(mat))){
    mat[i,j] <- mean(jaccard_index$jaccard_index[which(jaccard_index$label1==rownames(mat)[i] & jaccard_index$label2==colnames(mat)[j])])
  }
}

mat[mat == 1] <- NA  
pheatmap(mat,cluster_rows = F,cluster_cols = F,fontsize = 15,display_numbers = T,)

combinations <- paste(rep(tissues,each=length(conditions)*length(antibodys)),rep(antibodys, each = length(conditions)), rep(conditions, length(antibodys)), sep = "_")
mat <- matrix(NA, nrow = length(combinations), ncol = length(combinations), dimnames = list(combinations, combinations))  
antibody_group <- c(rep(paste(rep(antibodys, each = length(conditions)), rep(conditions, length(antibodys)), sep = "_"),times=length(tissues)))
tissue_group <- c(rep(tissues,each=6))
jaccard_index$tissue_label1 <- paste0(jaccard_index$tissue,"_",jaccard_index$antibody1,"_",jaccard_index$antibody1_condition)
jaccard_index$tissue_label2 <- paste0(jaccard_index$tissue,"_",jaccard_index$antibody2,"_",jaccard_index$antibody2_condition)

for(i in c(1:nrow(mat))){
  for(j in c(1:ncol(mat))){
    mat[i,j] <- mean(jaccard_index$jaccard_index[which(jaccard_index$tissue_label1==rownames(mat)[i] & jaccard_index$tissue_label2==colnames(mat)[j])])
    }
}

mat[mat == 1] <- NA  
annotation<-data.frame(antibody=antibody_group,tissue=tissue_group)
rownames(annotation)<-paste0(annotation$tissue,"_",annotation$antibody)

antibody_newCols <- colorRampPalette(hcl.colors(length(unique(annotation$antibody))))
antibody_annoCol <- antibody_newCols(length(unique(annotation$antibody)))
names(antibody_annoCol) <- unique(annotation$antibody)
tissue_newCols <- colorRampPalette(hcl.colors(length(unique(annotation$tissue)),palette = "Dynamic"))
tissue_annoCol <- tissue_newCols(length(unique(annotation$tissue)))
names(tissue_annoCol) <- unique(annotation$tissue)
ann_colors <- list(tissue=tissue_annoCol,antibody=antibody_annoCol)
pheatmap(mat,cluster_rows = F,cluster_cols = F,annotation_row = annotation,annotation_col = annotation,fontsize = 10,
         show_rownames = F,show_colnames = F,na_col = "white",border_color = "#dbe2ef", gaps_row =seq(6, 48, by = 6),gaps_col =seq(6, 48, by = 6),annotation_colors = ann_colors)


# jaccard_index$label <- factor(jaccard_index$label,levels=c(paste0(antibodys,"_",conditions)))
ggplot(jaccard_index[which(jaccard_index$antibody1_condition=="Up"),],mapping = aes(x=antibody1,y=jaccard_index,fill = label))+
  geom_boxplot()+theme_bw()+xlab("")+ylab("jaccard index")+ggtitle(paste0("Increase"))+
  theme(text = element_text(size = 18))
ggplot(jaccard_index[which(jaccard_index$antibody1_condition=="Down"),],mapping = aes(x=antibody1,y=jaccard_index,fill = label))+
  geom_boxplot()+theme_bw()+xlab("")+ylab("jaccard index")+ggtitle(paste0("Decrease"))+
  theme(text = element_text(size = 18))

jaccard_index<-data.frame(tissue = character(), 
                          jaccard_index = numeric(),
                          antibody1 = character(),
                          antibody1_condition=character(),
                          antibody2 = character(),
                          antibody2_condition=character(),
                          stringsAsFactors = FALSE)

conditions<-c("Up","Down")
antibodys <- c("H3K27ac","H3K4me3","H3K4me1")
for (i in c(1:length(tissues))){
  tissue<-tissues[i]
  for(j in c(1:length(antibodys))){
    antibody1 <- antibodys[j]
    for(k in c(1:length(antibodys))){
      antibody2 <- antibodys[k]
      df1 <- read.csv(paste0("data/samples/",tissue,"/",antibody1,"/",antibody1,"_",bin_size(antibody1),"_bins_diff_after_remove_batch_effect.csv")) 
      df2 <- read.csv(paste0("data/samples/",tissue,"/",antibody2,"/",antibody2,"_",bin_size(antibody2),"_bins_diff_after_remove_batch_effect.csv")) 
      for (condition1 in conditions){
        peaks1 <- df1$Geneid[which(df1$Significant_bar==condition1)]
        for (condition2 in conditions){
          peaks2 <- df2$Geneid[which(df2$Significant_bar==condition2)]
          t_jaccard_index <- data.frame(tissue = tissue, 
                                          jaccard_index = jaccard(peaks1,peaks2),
                                          antibody1 = antibody1,
                                          antibody1_condition=condition1,
                                          antibody2 = antibody2,
                                          antibody2_condition=condition2,
                                          stringsAsFactors = FALSE)
          jaccard_index <- rbind(jaccard_index,t_jaccard_index)
        }
      }
    }
  }
}

# jaccard_index$label <- paste0(jaccard_index$antibody2,"_",jaccard_index$antibody2_condition)
# jaccard_index$label <- factor(jaccard_index$label,levels=c(paste0(antibodys,"_",conditions)))
# ggplot(jaccard_index[which(jaccard_index$antibody1_condition=="Up"),],mapping = aes(x=antibody1,y=jaccard_index,fill = label))+
#   geom_boxplot()+theme_bw()+xlab("")+ylab("jaccard index")+ggtitle(paste0("Increase"))+
#   theme(text = element_text(size = 18))
# ggplot(jaccard_index[which(jaccard_index$antibody1_condition=="Down"),],mapping = aes(x=antibody1,y=jaccard_index,fill = label))+
#   geom_boxplot()+theme_bw()+xlab("")+ylab("jaccard index")+ggtitle(paste0("Decrease"))+
#   theme(text = element_text(size = 18))
# 

combinations <- paste(rep(antibodys, each = length(conditions)), rep(conditions, length(antibodys)), sep = "_")
mat <- matrix(0, nrow = length(combinations), ncol = length(combinations), dimnames = list(combinations, combinations))  
jaccard_index$label2 <- paste0(jaccard_index$antibody2,"_",jaccard_index$antibody2_condition)
jaccard_index$label1 <- paste0(jaccard_index$antibody1,"_",jaccard_index$antibody1_condition)
# jaccard_index[is.na(jaccard_index)] <- 0


for(i in c(1:nrow(mat))){
  for(j in c(1:ncol(mat))){
    mat[i,j] <- mean(jaccard_index$jaccard_index[which(jaccard_index$label1==rownames(mat)[i] & jaccard_index$label2==colnames(mat)[j])],na.rm = TRUE)
  }
}

mat[mat == 1] <- NA  
pheatmap(mat,cluster_rows = F,cluster_cols = F,fontsize = 15,display_numbers = T,)

combinations <- paste(rep(tissues,each=length(conditions)*length(antibodys)),rep(antibodys, each = length(conditions)), rep(conditions, length(antibodys)), sep = "_")
mat <- matrix(NA, nrow = length(combinations), ncol = length(combinations), dimnames = list(combinations, combinations))  
antibody_group <- c(rep(paste(rep(antibodys, each = length(conditions)), rep(conditions, length(antibodys)), sep = "_"),times=length(tissues)))
tissue_group <- c(rep(tissues,each=6))
jaccard_index$tissue_label1 <- paste0(jaccard_index$tissue,"_",jaccard_index$antibody1,"_",jaccard_index$antibody1_condition)
jaccard_index$tissue_label2 <- paste0(jaccard_index$tissue,"_",jaccard_index$antibody2,"_",jaccard_index$antibody2_condition)

for(i in c(1:nrow(mat))){
  for(j in c(1:ncol(mat))){
    mat[i,j] <- mean(jaccard_index$jaccard_index[which(jaccard_index$tissue_label1==rownames(mat)[i] & jaccard_index$tissue_label2==colnames(mat)[j])])
  }
}

mat[mat == 1] <- NA  
annotation<-data.frame(antibody=antibody_group,tissue=tissue_group)
rownames(annotation)<-paste0(annotation$tissue,"_",annotation$antibody)

antibody_newCols <- colorRampPalette(hcl.colors(length(unique(annotation$antibody))))
antibody_annoCol <- antibody_newCols(length(unique(annotation$antibody)))
names(antibody_annoCol) <- unique(annotation$antibody)
tissue_newCols <- colorRampPalette(hcl.colors(length(unique(annotation$tissue)),palette = "Dynamic"))
tissue_annoCol <- tissue_newCols(length(unique(annotation$tissue)))
names(tissue_annoCol) <- unique(annotation$tissue)
ann_colors <- list(tissue=tissue_annoCol,antibody=antibody_annoCol)
pheatmap(mat,cluster_rows = F,cluster_cols = F,annotation_row = annotation,annotation_col = annotation,fontsize = 10,
         show_rownames = F,show_colnames = F,na_col = "white",border_color = "#dbe2ef", gaps_row =seq(6, 48, by = 6),gaps_col =seq(6, 48, by = 6),annotation_colors = ann_colors)
