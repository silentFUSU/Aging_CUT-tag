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

all_10kb_bins <- read.delim("~/ref_data/mm10_10kb_bins.bed",header = F)
all_10kb_bins <- all_10kb_bins$V4
jaccard_index_random<-data.frame(tissue = character(), 
                                 jaccard_index = numeric(),
                                 antibody1 = character(),
                                 antibody1_condition=character(),
                                 antibody2 = character(),
                                 antibody2_condition=character(),
                                 stringsAsFactors = FALSE)
conditions<-c("Up","Down")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","Hip","cecum")
random_time <- 10000
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
          t_jaccard_index=0
          for(random in c(1:random_time)){
            peaks1 <- sample(all_10kb_bins,length(peaks1))
            peaks2 <- sample(all_10kb_bins,length(peaks2))
            t_jaccard_index=t_jaccard_index+jaccard(peaks1,peaks2)
          }
          t_jaccard_index=t_jaccard_index/random_time
          t_jaccard_index <- data.frame(tissue = tissue, 
                                        jaccard_index = jaccard(peaks1,peaks2),
                                        antibody1 = antibody1,
                                        antibody1_condition=condition1,
                                        antibody2 = antibody2,
                                        antibody2_condition=condition2,
                                        stringsAsFactors = FALSE)
          jaccard_index_random <- rbind(jaccard_index_random,t_jaccard_index)
        }
      }
    }
  }
}
write.csv(jaccard_index_random,"result/all/jaccard_index_random.csv",row.names = F)