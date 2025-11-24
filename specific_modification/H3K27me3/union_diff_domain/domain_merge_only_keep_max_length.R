rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
options(scipen = 999) 
library(tidyverse)  
library(dplyr) 
library(data.table)

domain_pool <- read.table("data/samples/all/H3K27me3/bed/H3K27me3_edd_domain_pool.bed")
domain_pool$length <- domain_pool$V3 - domain_pool$V2 +1
domain_pool <- domain_pool[order(domain_pool$length,decreasing = T),]
domain_pool$label <- paste(domain_pool$V1,domain_pool$V2,domain_pool$V3,sep = "-")

domain_summary <- data.frame()
while(nrow(domain_pool) > 0){
  t_domain_pool <- domain_pool
  t_domain_pool <- t_domain_pool[order(t_domain_pool$length,decreasing = T),]
  
  domain <- t_domain_pool[1,]
  domain <- as.data.table(domain)
  setDT(domain)
  setkey(domain,V1,V2,V3)
  
  t_domain_pool <- as.data.table(t_domain_pool)
  setDT(t_domain_pool)
  setkey(t_domain_pool,V1,V2,V3)
  
  overlaps <- foverlaps(t_domain_pool,domain, type = "any", nomatch = 0L)  
  domain_pool <- domain_pool[which(! domain_pool$label %in% overlaps$i.label),]
  domain_summary <- rbind(domain_summary,domain)
}

domain_summary$V1 <- factor(domain_summary$V1,levels=paste0("chr",c(1:19,"X","Y")))
domain_summary$V2 <- as.numeric(domain_summary$V2)
domain_summary$V3 <- as.numeric(domain_summary$V3)
domain_summary <- domain_summary %>%
  arrange(V1, V2, V3)
write.table(domain_summary[,c(1:3)],"data/samples/all/H3K27me3/bed/H3K27me3_edd_domain_keep_max.bed",append = F,quote = F,sep = "\t",row.names = F,col.names = F)
