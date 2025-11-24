rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
tissue <- "BAT"
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
    }else if(tissue_label == "Mammarygland"){
      tissue_label <- "Mammary Gland"
    }
  }
  return(tissue_label)
} 
H3K27me3_change_in_H3K9me3_decrease_region <- function(tissue){
  # H3K9me3 <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"))
  # H3K9me3_decrease_region <- H3K9me3[which(H3K9me3$Significant_bar=="Down"),]
  # write.table(H3K9me3_decrease_region[,c("Chr","Start","End","Geneid")], 
  #             paste0("data/samples/",tissue,"/H3K9me3/bed/H3K9me3_10kb_bins_diff_after_remove_batch_effect_down.bed"), 
  #             sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  dir.create(paste0("result/",tissue,"/H3K27me3_H3K9me3_relationship/matrix"))
  dir.create(paste0("result/",tissue,"/H3K27me3_H3K9me3_relationship/plot"))
  search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue_label_change(tissue) & search_table$antibody== "H3K27me3"),]
  young <-search_table$sample_name[which(search_table$age=="3m")]
  old <- search_table$sample_name[which(search_table$age=="24m")]
  bw_27me3_rep1 <- paste0("data/samples/",tissue,"/H3K27me3/bw/",young[1],"*.nodup.bw"," ","data/samples/",tissue,"/H3K27me3/bw/",old[1],"*.nodup.bw")
  bw_27me3_rep2 <- paste0("data/samples/",tissue,"/H3K27me3/bw/",young[2],"*.nodup.bw"," ","data/samples/",tissue,"/H3K27me3/bw/",old[2],"*.nodup.bw")
  search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
  search_table <- search_table[which(search_table$tissue == tissue_label_change(tissue) & search_table$antibody== "H3K9me3"),]
  young <-search_table$sample_name[which(search_table$age=="3m")]
  old <- search_table$sample_name[which(search_table$age=="24m")]
  bw_9me3_rep1 <- paste0("data/samples/",tissue,"/H3K9me3/bw/",young[1],"*.nodup.bw"," ","data/samples/",tissue,"/H3K9me3/bw/",old[1],"*.nodup.bw")
  bw_9me3_rep2 <- paste0("data/samples/",tissue,"/H3K9me3/bw/",young[2],"*.nodup.bw"," ","data/samples/",tissue,"/H3K9me3/bw/",old[2],"*.nodup.bw")
  script_compute_matrix <- paste0("cd /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/ \n",
  "source /storage/zhangyanxiaoLab/suzhuojie/miniconda3/etc/profile.d/conda.sh \n",
  "conda activate snakemake\n",
  "computeMatrix scale-regions -S ",bw_27me3_rep1," -R data/samples/",tissue,"/H3K9me3/bed/H3K9me3_10kb_bins_diff_after_remove_batch_effect_down.bed ",
  "--beforeRegionStartLength 10000 --startLabel start --endLabel end ",
  "--regionBodyLength 10000 ",
  "--afterRegionStartLength 10000 ",
  "--numberOfProcessors 20 ",
  "--skipZeros -o result/",tissue,"/H3K27me3_H3K9me3_relationship/matrix/H3K27me3_in_H3K9me3_decreae_rep1.mat.gz & \n",
  "computeMatrix scale-regions -S ",bw_27me3_rep2," -R data/samples/",tissue,"/H3K9me3/bed/H3K9me3_10kb_bins_diff_after_remove_batch_effect_down.bed ",
  "--beforeRegionStartLength 10000 --startLabel start --endLabel end ",
  "--regionBodyLength 10000 ",
  "--afterRegionStartLength 10000 ",
  "--numberOfProcessors 20 ",
  "--skipZeros -o result/",tissue,"/H3K27me3_H3K9me3_relationship/matrix/H3K27me3_in_H3K9me3_decreae_rep2.mat.gz & \n",
  "wait \n",
  "echo compute matrix done"
  )
  system(script_compute_matrix)
  script_plot <- paste0("cd /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/ \n",
  "source /storage/zhangyanxiaoLab/suzhuojie/miniconda3/etc/profile.d/conda.sh \n",
  "conda activate snakemake\n",
  "plotHeatmap -m result/",tissue,"/H3K27me3_H3K9me3_relationship/matrix/H3K27me3_in_H3K9me3_decreae_rep1.mat.gz ",
  "-o result/",tissue,"/H3K27me3_H3K9me3_relationship/plot/H3K27me3_in_H3K9me3_decrease_rep1.pdf ",
  "--colorMap viridis --startLabel start --endLabel end --perGroup --samplesLabel Young Old ",
  "--legendLocation none --plotTitle ",gsub(" ","_",tissue_label_change(tissue)),"_H3K27me3_in_H3K9me3_decrease & \n",
  "plotHeatmap -m result/",tissue,"/H3K27me3_H3K9me3_relationship/matrix/H3K27me3_in_H3K9me3_decreae_rep2.mat.gz ",
  "-o result/",tissue,"/H3K27me3_H3K9me3_relationship/plot/H3K27me3_in_H3K9me3_decrease_rep2.pdf ",
  "--colorMap viridis --startLabel start --endLabel end --perGroup --samplesLabel Young Old ",
  "--legendLocation none --plotTitle ",gsub(" ","_",tissue_label_change(tissue)),"_H3K27me3_in_H3K9me3_decrease & \n",
  "wait \n",
  "echo plot done"
  )
  system(script_plot)
}

tissue <- "mammarygland"
H3K27me3_change_in_H3K9me3_decrease_region(tissue)
