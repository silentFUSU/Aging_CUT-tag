rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(corrplot)  
tissue <- "brain"
antibodys <- c("H3K27ac","H3K4me3","H3K4me1","ATAC")

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
    }
  }
  return(tissue_label)
} 

dir.create("result/all/correlation_between_epi_markers/")
correlation_between_active_markers <- function(tissue,antibodys){
  for(i in c(1:length(antibodys))){
    antibody <- antibodys[i]
    if(antibody == "ATAC"){
      t_df <- read.csv(paste0("data/samples/ATAC/",tissue,"/",antibody,"/",antibody,"_1kb_bins_diff.csv"))
    }else{
      t_df <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_1kb_bins_diff.csv"))
    }
    t_df <- t_df[,c("Geneid","LogFC.old.young")]  
    colnames(t_df) <- c("Geneid",paste0(antibody
    ))
    if(i == 1){
      df <- t_df
    }else{
      df <- merge(df,t_df,by="Geneid")
    }
  }
  cor_matrix <- cor(df[,2:5])
  p_matrix <- cor.mtest(df[,2:5])$p
  png(paste0("result/all/correlation_between_epi_markers/",tissue,"_active_markers_correlation.png"),width = 2400, height = 1800, res = 300,type="cairo")  
  corrplot(cor_matrix, method = "color",
           tl.col="black",addgrid.col="white",p.mat=p_matrix, 
           insig = "label_sig",sig.level=0.05,type="lower",
           diag = FALSE)  
  mtext(tissue_label_change(tissue), side = 3, line = 2, adj=0.7, cex = 1.5)  
  dev.off()
}
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder","CB","jejunum","uterus","ovary")
tissues <- c("BAT","mammarygland")
for(tissue in tissues){
  correlation_between_active_markers(tissue,antibodys)
}
