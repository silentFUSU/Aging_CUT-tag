rm(list=ls()) 
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(rtracklayer)
library(gridExtra)
library(grid)  
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggsignif)
library(edgeR)
options(bitmapType="cairo")  
gtf <- import("~/ref_data/TE_reference/mm10_rmsk_TE.gtf", format = "gtf")
family_data <- as.data.frame(gtf[,c("gene_id","transcript_id","family_id","class_id")])
family_data <- as.data.table(family_data)
setDT(family_data)
setkey(family_data,seqnames,start,end)

regions <- read.csv(paste0("data/samples/all/H3K9me3/recursion_peaks_diff_table/kmeans_annotation_add_stable.csv"),row.names = 1)
regions <- regions %>%
  separate(label, into = c("chr", "range"), sep = ":") %>%
  separate(range, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end))

regions <- as.data.table(regions)
setDT(regions)
setkey(regions,chr,start,end)
overlaps <- as.data.frame(foverlaps(family_data, regions, type = "any", nomatch = 0L)) 
overlaps$cover_length <-pmin(overlaps$end,overlaps$i.end)- pmax(overlaps$i.start,overlaps$start)
overlaps <- as.data.frame(overlaps[,c("cluster","gene_id","transcript_id","family_id","class_id","cover_length")])
overlaps <- overlaps %>%
  group_by(gene_id, transcript_id, family_id, class_id) %>%
  slice_max(order_by = cover_length, n = 1, with_ties = FALSE) %>%
  ungroup()

overlaps$label <- paste(overlaps$transcript_id,overlaps$gene_id,overlaps$family_id,overlaps$class_id,sep=":")

TE_info <- as.data.frame(gtf[,c("gene_id","transcript_id","family_id","class_id")])
TE_info$label <- paste(TE_info$transcript_id,TE_info$gene_id,TE_info$family_id,TE_info$class_id,sep=":")
TE_info <- TE_info[,c("label","gene_id")]

tissues <- c("aorta","BAT","bladder","bonemarrow","brain","CB","cecum","colon","heart","Hip","ileum","jejunum","kidney","liver",
             "lung","muscle","ovary","pancreas","skin","spleen","stomach","testis","thymus","tongue","uterus","mammarygland","iWAT")
for(tissue in tissues){
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  search_table <- search_table[which(search_table$tissue_label==tissue),]
  tab <- read.table(paste0("data/samples/RNA/",tissue,"/TElocal/combined.cntTable"),header = T,row.names = 1)  
  counts <- tab
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SRR[0-9]+|HM[0-9]+).*"
  colnames(counts) <- gsub(pattern, "\\1", colnames(counts))  
  counts <- counts[,search_table$sample_name]
  counts$condition <- "other"
  counts$condition[which(rownames(counts) %in% overlaps$label[which(overlaps$cluster=="kmeans1")])] <- "kmeans1"
  counts$condition[which(rownames(counts) %in% overlaps$label[which(overlaps$cluster=="kmeans2")])] <- "kmeans2"
  counts$condition[which(rownames(counts) %in% overlaps$label[which(overlaps$cluster=="kmeans3")])] <- "kmeans3"
  counts$condition[which(rownames(counts) %in% overlaps$label[which(overlaps$cluster=="kmeans4")])] <- "kmeans4"
  counts$condition[which(rownames(counts) %in% overlaps$label[which(overlaps$cluster=="Stable")])] <- "Stable"
  
  genes <- rownames(counts)[grep("^ENSMUSG", rownames(counts))]
  genes_info <- data.frame(label=genes,gene_id=genes)
  genes_TE_info <- rbind(genes_info,TE_info)
  genes_TE_info$label <- factor(genes_TE_info$label,levels=rownames(counts))
  genes_TE_info <- genes_TE_info[order(genes_TE_info$label),]
  
  counts$gene_id <- genes_TE_info$gene_id
  
  counts <- counts %>%
    group_by(condition,gene_id) %>%
    summarise(
      across(.cols = 1:nrow(search_table), sum, .names = "{.col}")
    )
  counts <- as.data.frame(counts)
  rownames(counts) <- paste0(counts$gene_id,"_",counts$condition)
  counts <- counts[,search_table$sample_name]
  search_table$sample_name <- factor(search_table$sample_name,colnames(counts))
  search_table <- search_table[order(search_table$sample_name),]
  age <- search_table$age
  age[which(age=="3m")] <- "young"
  age[which(age=="24m")] <- "old"
  mouse_ID <- search_table$mouse_ID
  colnames(counts) <- paste0(colnames(counts),"-",mouse_ID,"-",age)
  
  y= DGEList(counts=counts,group=age)
  keep = which(rowSums(cpm(y)>0)>=2)
  y = y[keep,]
  y$samples$group <- factor(y$samples$group, levels=c("young","old"))
  design <- model.matrix(~group, y$samples)
  y <- calcNormFactors(y)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, coef = 2)
  out = cbind(cpm(y),lrt$table, "fdr"=p.adjust(lrt$table$PValue,method="BH"))
  out$Significant <- ifelse(out$fdr< 0.05 & abs(out$logFC) >= 0, 
                            ifelse(out$logFC > 0, "Up", "Down"), "Stable")
  
  out <- out[-grep("^ENSMUSG", rownames(out)),]
  write.csv(out,paste0("data/samples/RNA/",tissue,"/diff_expression_TE_local_H3K9me3_kmeans_change_filter_bar.csv"))
}
