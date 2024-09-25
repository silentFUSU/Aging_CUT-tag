rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
args <- commandArgs(trailingOnly = TRUE)  
if (length(args) < 1) {  
  stop("No tissue argument provided")  
}  
tissue <- args[1]  
print(paste("Tissue is:", tissue))  

make_tss_ref <- function(tissue){
  tss_ref <- read.table("~/ref_data/for_normal_mapping/mm10/mm10_refseq_TSS.bed")
  diff_gene <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene.csv"))
  increase_gene <- diff_gene$X[which(diff_gene$Significant=="Up")]
  decrease_gene <- diff_gene$X[which(diff_gene$Significant=="Down")]
  tss_increase <- tss_ref[tss_ref$V6 %in% increase_gene,]
  tss_decrease <- tss_ref[tss_ref$V6 %in% decrease_gene,]
  dir.create(paste0("data/samples/RNA/",tissue,"/bed"))
  write.table(tss_increase,paste0("data/samples/RNA/",tissue,"/bed/increase_gene_TSS.bed"),col.names = F,row.names = F,quote = F,sep = "\t")
  write.table(tss_decrease,paste0("data/samples/RNA/",tissue,"/bed/decrease_gene_TSS.bed"),col.names = F,row.names = F,quote = F,sep = "\t")
}
make_tss_ref(tissue)