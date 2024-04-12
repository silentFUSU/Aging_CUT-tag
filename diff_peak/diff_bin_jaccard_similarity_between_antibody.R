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
library(pheatmap)
library(reshape2)
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

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# jaccard_index<-data.frame(tissue = character(),
#                           jaccard_index = numeric(),
#                           antibody1 = character(),
#                           antibody1_condition=character(),
#                           antibody2 = character(),
#                           antibody2_condition=character(),
#                           stringsAsFactors = FALSE)
jaccard_index <-read.csv("result/all/jaccard_index.csv")
conditions<-c("Up","Down")
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")

tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","Hip","cecum","bonemarrow","ileum","heart","thymus")
# tissues <- c("bonemarrow")
# tissues <- c("ileum")
tissues <- c("heart","thymus")
pb <- txtProgressBar(min = 0, length(tissues)*length(antibodys)*length(antibodys), style = 3)
counter=0
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
      counter <- counter+1
      setTxtProgressBar(pb, counter)
    }
  }
}
# write.csv(jaccard_index,"result/all/jaccard_index.csv",row.names = F)


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
pl <- pheatmap(mat,cluster_rows = F,cluster_cols = F,fontsize = 15,display_numbers = T, breaks = seq(0, 0.09, length.out = 101))
save_pheatmap_pdf(pl,paste0("result/all/jaccard_index/all_heatmap_bin_in_peak.pdf"),height=7,width = 8)

jaccard_index_random <-read.csv("result/all/jaccard_index_random.csv")

mat_random <- matrix(0, nrow = length(combinations), ncol = length(combinations), dimnames = list(combinations, combinations))  
jaccard_index_random$label2 <- paste0(jaccard_index_random$antibody2,"_",jaccard_index_random$antibody2_condition)
jaccard_index_random$label1 <- paste0(jaccard_index_random$antibody1,"_",jaccard_index_random$antibody1_condition)
jaccard_index_random$jaccard_index[which(jaccard_index_random$antibody1==jaccard_index_random$antibody2)]<-0
jaccard_index_random$jaccard_index[which(jaccard_index_random$label1==jaccard_index_random$label2)]<-1
for(i in c(1:nrow(mat_random))){
  for(j in c(1:ncol(mat_random))){
    mat_random[i,j] <- mean(jaccard_index_random$jaccard_index[which(jaccard_index_random$label1==rownames(mat_random)[i] & jaccard_index_random$label2==colnames(mat_random)[j])],na.rm = TRUE)
  }
}
mat_random[mat_random == 1] <- NA  
for(i in c(1:(nrow(mat_random)-1))){
  for(j in c((i+1):ncol(mat_random))){
    mat_random[i,j] <- mat_random[j,i]
  }
}

colors <- colorRampPalette(c("#07689f", "#ffde7d", "#e84545"))(100) 

pl <- pheatmap(mat_random,cluster_rows = F,cluster_cols = F,fontsize = 15,display_numbers = T, breaks = seq(0, 0.09, length.out = 101))
save_pheatmap_pdf(pl,paste0("result/all/jaccard_index/all_heatmap_random.pdf"),height=7,width = 8)


combinations <- paste(rep(antibodys, each = length(conditions)), rep(conditions, length(antibodys)), sep = "_")
mat_compare <- matrix(NA, nrow = length(combinations), ncol = length(combinations), dimnames = list(combinations, combinations))
for(i in c(1:nrow(mat_compare))){
  for(j in c(1:i)){
    if(rownames(mat_compare)[i] != colnames(mat_compare)[j]){
      df <- merge(jaccard_index[which(jaccard_index$label1==rownames(mat_compare)[i] &jaccard_index$label2==colnames(mat_compare)[j]),c("tissue","jaccard_index")],
                  jaccard_index_random[which(jaccard_index_random$label1==rownames(mat_compare)[i]&jaccard_index_random$label2==colnames(mat_compare)[j]),c("tissue","jaccard_index")],
                  by="tissue")
      colnames(df)[2:3] <- c("reality","random")
      sig <- wilcox.test(df$reality,df$random, paired = TRUE, alternative = "greater")
      if(sig$p.value < 0.05){
        mat_compare[i,j] <- -log10(sig$p.value)
      }
    }
  }
}
for(i in c(1:(nrow(mat_compare)-1))){
  for(j in c((i+1):ncol(mat_compare))){
    mat_compare[i,j] <- mat_compare[j,i]
  }
}

pl <- pheatmap(mat_compare,cluster_rows = F,cluster_cols = F,fontsize = 15,breaks = seq(-log10(0.05), 3, length.out = 101))
save_pheatmap_pdf(pl,paste0("result/all/jaccard_index/all_heatmap_compare.pdf"),height=7,width = 8)

##########################################################
jaccard_index <-read.csv("result/all/jaccard_index.csv") 
jaccard_index_random <-read.csv("result/all/jaccard_index_random.csv")
jaccard_index <- read.csv("result/all/jaccard_index_bin_in_peak.csv")
jaccard_index_random$label2 <- paste0(jaccard_index_random$antibody2,"_",jaccard_index_random$antibody2_condition)
jaccard_index_random$label1 <- paste0(jaccard_index_random$antibody1,"_",jaccard_index_random$antibody1_condition)
jaccard_index$label2 <- paste0(jaccard_index$antibody2,"_",jaccard_index$antibody2_condition)
jaccard_index$label1 <- paste0(jaccard_index$antibody1,"_",jaccard_index$antibody1_condition)

H3K27me3_up_H3K9me3_down <- merge(jaccard_index[which(jaccard_index$label1=="H3K27me3_Up"&jaccard_index$label2=="H3K9me3_Down"),c("tissue","jaccard_index")],
                                  jaccard_index_random[which(jaccard_index_random$label1=="H3K27me3_Up"&jaccard_index_random$label2=="H3K9me3_Down"),c("tissue","jaccard_index")],by="tissue")
wilcox.test(H3K27me3_up_H3K9me3_down$jaccard_index.x,H3K27me3_up_H3K9me3_down$jaccard_index.y, paired = TRUE, alternative = "greater")

colnames(H3K27me3_up_H3K9me3_down)[2:3] <- c("reality","random")
H3K27me3_up_H3K9me3_down_long<- melt(H3K27me3_up_H3K9me3_down, id.vars = c("tissue"), measure.vars = c("reality","random"), variable.name = "condition", value.name = "jaccard_index")  
ggplot()+
  geom_boxplot(H3K27me3_up_H3K9me3_down_long,mapping=aes(x=condition,y=jaccard_index,fill =condition))+
  geom_jitter(H3K27me3_up_H3K9me3_down_long,mapping=aes(x=condition,y=jaccard_index,color =tissue),position = position_jitter(width = 0.3), size = 3, alpha = 0.7)+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("Peaks")+
  xlab("")+labs(fill = "", color = "")

H3K27me3_down_H3K9me3_up <- merge(jaccard_index[which(jaccard_index$label1=="H3K27me3_Down"&jaccard_index$label2=="H3K9me3_Up"),c("tissue","jaccard_index")],
                                  jaccard_index_random[which(jaccard_index_random$label1=="H3K27me3_Down"&jaccard_index_random$label2=="H3K9me3_Up"),c("tissue","jaccard_index")],by="tissue")
wilcox.test(H3K27me3_down_H3K9me3_up$jaccard_index.x,H3K27me3_down_H3K9me3_up$jaccard_index.y, paired = TRUE, alternative = "greater")

colnames(H3K27me3_down_H3K9me3_up)[2:3] <- c("reality","random")
H3K27me3_down_H3K9me3_up_long<- melt(H3K27me3_down_H3K9me3_up, id.vars = c("tissue"), measure.vars = c("reality","random"), variable.name = "condition", value.name = "jaccard_index")  
ggplot()+
  geom_boxplot(H3K27me3_down_H3K9me3_up_long,mapping=aes(x=condition,y=jaccard_index,fill =condition))+
  geom_jitter(H3K27me3_down_H3K9me3_up_long,mapping=aes(x=condition,y=jaccard_index,color =tissue),position = position_jitter(width = 0.3), size = 3, alpha = 0.7)+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("Peaks")+
  xlab("")+labs(fill = "", color = "")
pb <- txtProgressBar(min = 0, length(tissues), style = 3)
for(k in c(1:length(tissues))){
  plot_jaccard_index <- jaccard_index[which(jaccard_index$tissue==tissues[k]),]
  mat_tissue<- matrix(0, nrow = length(combinations), ncol = length(combinations), dimnames = list(combinations, combinations)) 
  for(i in c(1:nrow(mat_tissue))){
    for(j in c(1:ncol(mat_tissue))){
      mat_tissue[i,j] <- mean(plot_jaccard_index$jaccard_index[which(plot_jaccard_index$label1==rownames(mat_tissue)[i] & plot_jaccard_index$label2==colnames(mat_tissue)[j])],na.rm = TRUE)
    }
  }
  mat_tissue[mat_tissue == 1] <- NA  
  pl <-pheatmap(mat_tissue,cluster_rows = F,cluster_cols = F,fontsize = 15,display_numbers = T, breaks = seq(0, 0.1, length.out = 101))
  save_pheatmap_pdf(pl,paste0("result/all/jaccard_index/",tissues[k],"_bin_in_peak_heatmap.pdf"),height=7,width = 8)
  setTxtProgressBar(pb, k)
}


H3K9me3_down_H3K27ac_up <- merge(jaccard_index[which(jaccard_index$label1=="H3K9me3_Down"&jaccard_index$label2=="H3K27ac_Up"),c("tissue","jaccard_index")],
                                  jaccard_index_random[which(jaccard_index_random$label1=="H3K9me3_Down"&jaccard_index_random$label2=="H3K27ac_Up"),c("tissue","jaccard_index")],by="tissue")
wilcox.test(H3K9me3_down_H3K27ac_up$jaccard_index.x,H3K9me3_down_H3K27ac_up$jaccard_index.y, paired = TRUE, alternative = "greater")

colnames(H3K9me3_down_H3K27ac_up)[2:3] <- c("reality","random")
H3K9me3_down_H3K27ac_up_long<- melt(H3K9me3_down_H3K27ac_up, id.vars = c("tissue"), measure.vars = c("reality","random"), variable.name = "condition", value.name = "jaccard_index")  
ggplot()+
  geom_boxplot(H3K9me3_down_H3K27ac_up_long,mapping=aes(x=condition,y=jaccard_index,fill =condition))+
  geom_jitter(H3K9me3_down_H3K27ac_up_long,mapping=aes(x=condition,y=jaccard_index,color =tissue),position = position_jitter(width = 0.3), size = 3, alpha = 0.7)+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("Peaks")+
  xlab("")+labs(fill = "", color = "")


H3K9me3_down_H3K4me1_up <- merge(jaccard_index[which(jaccard_index$label1=="H3K9me3_Down"&jaccard_index$label2=="H3K4me1_Up"),c("tissue","jaccard_index")],
                                 jaccard_index_random[which(jaccard_index_random$label1=="H3K9me3_Down"&jaccard_index_random$label2=="H3K4me1_Up"),c("tissue","jaccard_index")],by="tissue")
wilcox.test(H3K9me3_down_H3K4me1_up$jaccard_index.x,H3K9me3_down_H3K4me1_up$jaccard_index.y, paired = TRUE, alternative = "greater")

colnames(H3K9me3_down_H3K4me1_up)[2:3] <- c("reality","random")
H3K9me3_down_H3K4me1_up_long<- melt(H3K9me3_down_H3K4me1_up, id.vars = c("tissue"), measure.vars = c("reality","random"), variable.name = "condition", value.name = "jaccard_index")  
ggplot()+
  geom_boxplot(H3K9me3_down_H3K4me1_up_long,mapping=aes(x=condition,y=jaccard_index,fill =condition))+
  geom_jitter(H3K9me3_down_H3K4me1_up_long,mapping=aes(x=condition,y=jaccard_index,color =tissue),position = position_jitter(width = 0.3), size = 3, alpha = 0.7)+
  # geom_jitter(shape=16,size=1,position = position_jitter(0.2))+
  scale_fill_brewer(palette="Set3")+
  theme_bw()+theme(text = element_text(size = 18))+ylab("Peaks")+
  xlab("")+labs(fill = "", color = "")




