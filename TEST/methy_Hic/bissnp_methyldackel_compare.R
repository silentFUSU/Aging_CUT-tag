rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
# bissnp <- read.delim("data/raw_data/20240725_methylHiC/bhmem/WJH_Mousebrain_C1_BSseq_101/vcf/WJH_Mousebrain_C1_BSseq_101.calmd.nodup.cpg.raw.sort.CG.6plus2.bed",skip=1,header=F)
bissnp_mincov0 <- read.delim("data/raw_data/20240725_methylHiC/bhmem/WJH_Mousebrain_C1_BSseq_101/vcf_minconv0/WJH_Mousebrain_C1_BSseq_101.calmd.nodup.cpg.raw.sort.CG.6plus2.bed",skip=1,header=F)
# methyldackel <- read.delim("data/raw_data/20240725_methylHiC/bhmem/WJH_Mousebrain_C1_BSseq_101/vcf/WJH_Mousebrain_C1_BSseq_101.calmd.nodup_CpG.bedGraph",skip=1,header=F)
methyldackel_30 <- read.delim("data/raw_data/20240725_methylHiC/bhmem/WJH_Mousebrain_C1_BSseq_101/bam/WJH_Mousebrain_C1_BSseq_101.calmd.nodup_CpG.bedGraph",skip=1,header=F)
head(bissnp_mincov0)
head(methyldackel_30)

methyldackel_30$depth <- methyldackel_30$V5+methyldackel_30$V6
bissnp_mincov0$depth <- bissnp_mincov0$V8

methyldackel_30$label <- paste0(methyldackel_30$V1,"-",methyldackel_30$V2,"-",methyldackel_30$V3)
bissnp_mincov0$label <- paste0(bissnp_mincov0$V1,"-",bissnp_mincov0$V2,"-",bissnp_mincov0$V3)
df <- merge(methyldackel_30[,c("depth","label")],bissnp_mincov0[,c("depth","label")],by="label")
head(df)
df$depth.x <- as.numeric(df$depth.x)
df$depth.y <- as.numeric(df$depth.y)

df_order <- df[order(-df[,2]),]

head(df_order)
df_same <- df[which(df$depth.x==df$depth.y),]
df_same <- df_same[order(-df_same$depth.x),]

head(df_same)
