rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(ChIPseeker)
library(EnsDb.Mmusculus.v79)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
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
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","Hip","cecum")
for (i in c(9:length(tissues))){
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

## mean value of all tissues
combinations <- paste(rep(antibodys, each = length(conditions)), rep(conditions, length(antibodys)), sep = "_")
mat <- matrix(0, nrow = length(combinations), ncol = length(combinations), dimnames = list(combinations, combinations))  
jaccard_index$label2 <- paste0(jaccard_index$antibody2,"_",jaccard_index$antibody2_condition)
jaccard_index$label1 <- paste0(jaccard_index$antibody1,"_",jaccard_index$antibody1_condition)
for(i in c(1:nrow(mat))){
  for(j in c(1:ncol(mat))){
    mat[i,j] <- mean(jaccard_index$jaccard_index[which(jaccard_index$label1==rownames(mat)[i] & jaccard_index$label2==colnames(mat)[j])],na.rm = TRUE)
  }
}
mat[mat == 1] <- NA  
colors <- colorRampPalette(c("#07689f", "#ffde7d", "#e84545"))(100) 
pheatmap(mat,cluster_rows = F,cluster_cols = F,fontsize = 15,display_numbers = T, breaks = seq(0, 0.09, length.out = 101))

# all_10kb_bins <- read.delim("~/ref_data/mm10_10kb_bins.bed",header = F)
# all_10kb_bins <- all_10kb_bins$V4
# jaccard_index_random<-data.frame(tissue = character(), 
#                           jaccard_index = numeric(),
#                           antibody1 = character(),
#                           antibody1_condition=character(),
#                           antibody2 = character(),
#                           antibody2_condition=character(),
#                           stringsAsFactors = FALSE)
# for (i in c(1:length(tissues))){
#   tissue<-tissues[i]
#   for(j in c(1:length(antibodys))){
#     antibody1 <- antibodys[j]
#     for(k in c(1:length(antibodys))){
#       antibody2 <- antibodys[k]
#       df1 <- read.csv(paste0("data/samples/",tissue,"/",antibody1,"/",antibody1,"_",bin_size(antibody1),"_bins_diff_after_remove_batch_effect.csv")) 
#       df2 <- read.csv(paste0("data/samples/",tissue,"/",antibody2,"/",antibody2,"_",bin_size(antibody2),"_bins_diff_after_remove_batch_effect.csv")) 
#       for (condition1 in conditions){
#         peaks1 <- df1$Geneid[which(df1$Significant_bar==condition1)]
#         for (condition2 in conditions){
#           peaks2 <- df2$Geneid[which(df2$Significant_bar==condition2)]
#           t_jaccard_index=0
#           for(random in c(1:100)){
#             peaks1 <- sample(all_10kb_bins,length(peaks1))
#             peaks2 <- sample(all_10kb_bins,length(peaks2))
#             t_jaccard_index=t_jaccard_index+jaccard(peaks1,peaks2)
#           }
#           t_jaccard_index=t_jaccard_index/100
#           t_jaccard_index <- data.frame(tissue = tissue, 
#                                         jaccard_index = jaccard(peaks1,peaks2),
#                                         antibody1 = antibody1,
#                                         antibody1_condition=condition1,
#                                         antibody2 = antibody2,
#                                         antibody2_condition=condition2,
#                                         stringsAsFactors = FALSE)
#           jaccard_index_random <- rbind(jaccard_index_random,t_jaccard_index)
#         }
#       }
#     }
#   }
# }
# write.csv(jaccard_index_random,"result/all/jaccard_index_random.csv",row.names = F)
jaccard_index_random$jaccard_index[which(jaccard_index_random$antibody1==jaccard_index_random$antibody2)]<-0
jaccard_index_random$label2 <- paste0(jaccard_index_random$antibody2,"_",jaccard_index_random$antibody2_condition)
jaccard_index_random$label1 <- paste0(jaccard_index_random$antibody1,"_",jaccard_index_random$antibody1_condition)
jaccard_index_random$jaccard_index[which(jaccard_index_random$label1==jaccard_index_random$label2)]<-NA

mat_random <- matrix(0, nrow = length(combinations), ncol = length(combinations), dimnames = list(combinations, combinations))  
for(i in c(1:nrow(mat_random))){
  for(j in c(1:ncol(mat_random))){
    mat_random[i,j] <- mean(jaccard_index_random$jaccard_index[which(jaccard_index_random$label1==rownames(mat_random)[i] & jaccard_index_random$label2==colnames(mat_random)[j])],na.rm = TRUE)
  }
}
pheatmap(mat_random,cluster_rows = F,cluster_cols = F,fontsize = 15,display_numbers = T, breaks = seq(0, 0.09, length.out = 101))

H3K27me3_up_H3K9me3_down <- merge(jaccard_index[which(jaccard_index$label1=="H3K27me3_Up"&jaccard_index$label2=="H3K9me3_Down"),c("tissue","jaccard_index")],jaccard_index_random[which(jaccard_index_random$label1=="H3K27me3_Up"&jaccard_index_random$label2=="H3K9me3_Down"),c("tissue","jaccard_index")],by="tissue")
t.test(H3K27me3_up_H3K9me3_down$jaccard_index.x,H3K27me3_up_H3K9me3_down$jaccard_index.y)
