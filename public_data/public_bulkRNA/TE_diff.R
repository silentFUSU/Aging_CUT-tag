rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(tidyr) 
library(dplyr)
library(stringr)
library(reshape2)
library(biomaRt)  
gtf <- "/storage/zhangyanxiaoLab/share/gtf/mm10.gencode.vM25.annotation.gtf"
gtf_lines <- readLines(gtf)
gene_id_list <- list()  
gene_name_list <- list()  
for (line in gtf_lines) {  
  if (grepl("^#", line)) next  # 跳过注释行  
  fields <- strsplit(line, "\t")[[1]]  
  if (fields[3] == "gene") {  # 仅处理类型为 gene 的行  
    attributes <- strsplit(fields[9], ";")[[1]]  
    gene_id <- sub('gene_id "([^"]+)".*', '\\1', attributes[grep('gene_id', attributes)])  
    gene_name <- sub('.*gene_name "([^"]+)".*', '\\1', attributes[grep('gene_name', attributes)])  
    if (length(gene_id) > 0 && length(gene_name) > 0) {  
      gene_id_list[[gene_id]] <- gene_name  
    }  
  }  
}  
gene_id_vector <- as.vector(names(gene_id_list))  
gene_name_vector <- as.vector(unlist(gene_id_list)) 
gene_id_map <- setNames(gene_name_vector, gene_id_vector) 

search_table <- read.csv("data/public_data/GSE132040/GSE132040_MACA_Bulk_metadata.csv")
colnames(search_table)[11] <- "sample"
search_table$source.name <- sub("_(\\w+)_\\d*|_(\\w+)$", "\\1",  search_table$source.name)  
search_table <- search_table[which(search_table$source.name=="Skin" & search_table$characteristics..age %in% c("3","24") & search_table$characteristics..sex=="m"),]

tab <- read.delim(paste0("data/public_data/GSE132040/te_counts.txt"),row.names = 1)
counts <- tab
pattern <- ".*te\\.(SRR[0-9]+).*"
colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
counts <- counts[,which(colnames(counts) %in% search_table$sample)]
search_table$sample <- factor(search_table$sample, levels=colnames(counts))
search_table <- search_table[order(search_table$sample),]

age <- search_table$characteristics..age
colnames(counts) <- paste0(colnames(counts),"-",age)
age[which(age=="3")] <- "young"
age[which(age=="24")] <- "old"

y= DGEList(counts=counts,group=age)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y$samples$group <- factor(y$samples$group, levels=c("young","old"))
design <- model.matrix(~group, y$samples)
y <- calcNormFactors(y)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 2)
tab<-tab[keep,]

out = cbind(cpm(y),lrt$table, "fdr"=p.adjust(lrt$table$PValue,method="BH"))

out$Significant <- ifelse(out$fdr< 0.05 & abs(out$logFC) >= 0, 
                          ifelse(out$logFC > 0, "Up", "Down"), "Stable")
out$Significant[which(!grepl("^ENSMUSG", rownames(out)) & out$Significant=="Down") ] <- "TE-Down"
out$Significant[which(!grepl("^ENSMUSG", rownames(out)) & out$Significant=="Up") ] <- "TE-Up"
colour <- setNames(c("grey","#3490de","#ea5455", "blue","red"),c("Stable","Down","Up","TE-Down","TE-Up"))
top_TE_Down <- out %>%  
  filter(Significant == "TE-Down") %>%  
  arrange(fdr) %>%  
  head(20)  
split_names <- strsplit(rownames(top_TE_Down), ":")  
split_df <- do.call(rbind, split_names) 
top_TE_Down <- cbind(top_TE_Down, split_df)
colnames(top_TE_Down)[ncol(top_TE_Down)-2] <- "gene"

top_TE_Up <- out %>%  
  filter(Significant == "TE-Up") %>%  
  arrange(fdr) %>%  
  head(20)  
split_names <- strsplit(rownames(top_TE_Up), ":")  
split_df <- do.call(rbind, split_names) 
top_TE_Up <- cbind(top_TE_Up, split_df)
colnames(top_TE_Up)[ncol(top_TE_Up)-2] <- "gene"
out <- out[which(!grepl("^ENSMUSG", rownames(out))), ]
p <-ggplot() +
  geom_point(data=out[which(out$Significant=="Stable"),], mapping=aes( logFC,  -log10(fdr),color = Significant), size=2)+
  geom_point(data=out[which(out$Significant=="Down"),], mapping=aes( logFC,  -log10(fdr),color = Significant), size=2) +  
  geom_point(data=out[which(out$Significant=="Up"),], mapping=aes( logFC,  -log10(fdr),color = Significant), size=2) +
  geom_point(data=out[which(out$Significant=="TE-Down"),], mapping=aes( logFC,  -log10(fdr),color = Significant), size=2) +
  geom_point(data=out[which(out$Significant=="TE-Up"),], mapping=aes( logFC,  -log10(fdr),color = Significant), size=2) +
  scale_color_manual(values = colour) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (fdr)") +
  theme_bw()+
  theme(text = element_text(size = 20))+
  ggtitle("TMS-Skin")
p <- p+geom_text_repel(data=top_TE_Down, mapping=aes(x=logFC, y=-log10(fdr), label=gene), vjust=-1, size=3) +  
  geom_text_repel(data=top_TE_Up, mapping=aes(x=logFC, y=-log10(fdr), label=gene), vjust=-1, size=3)  
