rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)  
library(parallel)  
args <- commandArgs(trailingOnly = TRUE)  
if (length(args) < 1) {  
  stop("No tissue argument provided")  
}  
tissue <- args[1]  
bin_size <- args[2]
print(paste("Tissue is:", tissue))  
print(paste("bin size is:", bin_size))  

dir.create(paste0("data/samples/WGBS/",tissue,"/compress2bin/"))
search_table <- read.csv("data/samples/all/WGBS_search_table.csv")


compress_to_bin <- function(tissue, bin_size){
  samples <- search_table$sample_name[which(search_table$tissue==tissue)]
  bin_summary <- data.frame(bin = character(),
                            total_V4 = numeric(),
                            total_V5 = numeric(),
                            percent = numeric(),
                            tissue = character(),
                            sample = character())
  for(sample in samples){
    df <- fread(paste0("data/samples/WGBS/",tissue,"/bdg/",sample,"_CpG.bdg"),sep = "\t")
    df <- df[which(df$V1 %in% paste0("chr",c(1:19,"X","Y"))),]
    df$V1 <- factor(df$V1,paste0("chr",c(1:19,"X","Y")))
    df$V3 <- df$V2
    df <- df[order(df$V1,df$V2),]
    
    ref <- fread(paste0("~/ref_data/mm10_", bin_size, "_bins.bed")) 
    ref$V1 <- factor(ref$V1,paste0("chr",c(1:19,"X","Y")))
    ref <- ref[order(ref$V1,ref$V2),]
    ref$V2 <- ref$V2+1
    ref[, `:=`(total_V4 = NA, total_V5 = NA)]  
    
    setDT(df)  
    setDT(ref)  
    setkey(df, V1, V2, V3) 
    setkey(ref, V1, V2, V3) 
    overlaps <- foverlaps(df, ref, type = "any", nomatch = 0L)  
    results <- overlaps[, .(total_V4 = sum(i.V4), total_V5 = sum(V5)), by = .(V4)]  
    results$percent <- results$total_V4/results$total_V5
    colnames(results)[1] <- "label"
    results$tissue <- tissue
    results$sample <- sample
    bin_summary <- rbind(bin_summary,results)
    # t_df <- df %>% filter(V2 > "27570001" & V2 <= "27580000" & V1 == "chr1") #check
  }
  write.csv(bin_summary,paste0("data/samples/WGBS/",tissue,"/compress2bin/",bin_size,"_bins_all_depth.csv"),row.names = F)
}
compress_to_bin(tissue,bin_size)