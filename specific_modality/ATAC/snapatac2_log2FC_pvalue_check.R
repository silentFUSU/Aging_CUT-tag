rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(bitmapType="cairo")  
library(grid)
extract_before_bracket <- function(s) {  
  parts <- strsplit(s, "\\(")[[1]]  
  return(parts[1])  
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
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
tissue <- "lung"
condition <- "up"
data_path <- "data/samples/ATAC/ATAC_peak_from_LMJ/motif_snapatac2_cisbp/motif_bg/"
motif <- read.csv(paste0(data_path,condition,"/enrichment_results_",tissue,"_",str_to_title(condition),"_sorted.bed.csv"))

to_plot <- motif[,c("log2.fold.change.","adjusted.p.value")]
to_plot$adjusted.p.value[which(to_plot$adjusted.p.value == 0)] <- 1e-20

to_plot$adjusted.p.value <- -log10(to_plot$adjusted.p.value)
to_plot$log2.fold.change.[which(to_plot$log2.fold.change.== -Inf)] <- 0

ggplot(to_plot[which(to_plot$log2.fold.change.>0),], aes(x=`log2.fold.change.`, y=`adjusted.p.value`)) + 
  geom_point(size=1) +theme_bw()+
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  theme(text = element_text(size = 20))


