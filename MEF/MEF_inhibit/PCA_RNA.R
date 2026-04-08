rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
tab1 = read.delim(paste0("data/samples/RNA/MEF/combined-chrM.counts"),skip=1)
tab2 = read.delim(paste0("data/samples/RNA/MEF_EZH2_inhibit/combined-chrM.counts"),skip=1)

tab <- merge(tab1,tab2[,c(1,7:10)],by="Geneid")

tab <- tab[!grepl("chrY", tab$Chr), ]
rownames(tab) <- tab$Geneid

tab <- tab[,-1]
colnames <- colnames(tab)[6:length(tab)]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|XM[0-9]+|DYQ[0-9]+|TK[0-9]+).*"
new_colnames <- gsub(pattern, "\\1", colnames)
colnames(tab)[6:length(tab)] <- new_colnames
sorted_index <- order(new_colnames)
order_colnames <- new_colnames[sorted_index] 
counts <- tab[,order_colnames]   
group <- read.csv("data/samples/RNA/sample_tissue_info.csv",sep = ',')
group <- group[which(group$SampleID %in% colnames(counts)),]
group <- group[order(group$SampleID),]
group <- group[which(group$TissueName%in% c("MEF","MEF_mid","MEF_EZH2_inhibit")),]
counts <- counts[,group$SampleID]
age <- group[which(group$SampleID %in% colnames(counts)),"Age"]
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
search_table <- search_table[which(search_table$tissue_label %in% c("MEF","MEF_mid_age","MEF_EZH2_inhibit")),]
# search_table <- search_table[-c(9,10),]
# search_table <- search_table[-c(1,2),]

y= DGEList(counts=counts,group=age)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]

logCPMs <- cpm(y, log = TRUE)
pca <- prcomp(t(logCPMs))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)

to_plot <- merge(to_plot,search_table,by="sample_name")
to_plot$rownames <- paste0(to_plot$sample_name,"-",to_plot$mouse_ID,"-",to_plot$age)
to_plot$age <- factor(to_plot$age,levels=c("3m","24m"))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

p <- ggplot(to_plot, aes(x=PC1, y=PC2, color=age)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))+
  geom_text_repel(  
    data = to_plot,  
    aes(x = PC1, y = PC2, label = rownames, color = age),  
    size = 5,  
    box.padding = unit(0.35, "lines"),  
    point.padding = unit(0.3, "lines")  
  ) 
p

logCPMs_corrected <- limma::removeBatchEffect(logCPMs, batch = c("batch1","batch1","batch1","batch1","batch1","batch1","batch2","batch2","batch2","batch2"))

pca <- prcomp(t(logCPMs_corrected))
to_plot <- data.frame(pca$x)
to_plot$sample_name <- rownames(to_plot)
to_plot <- merge(to_plot,search_table,by="sample_name")
to_plot$rownames <- paste0(to_plot$sample_name,"-",to_plot$mouse_ID,"-",to_plot$age)
to_plot$age <- factor(to_plot$age,levels=c("3m","24m"))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))


to_plot$age <- c("p2","p2","p6","p6","p10","p10","DMSO","DMSO","GSK126","GSK126")
to_plot$condition <- c("RS","RS","RS","RS","RS","RS","EZH2i","EZH2i","EZH2i","EZH2i")
p <- ggplot(to_plot, aes(x=PC1, y=PC2, color=age, shape=condition)) + 
  geom_point(size=5) +theme_bw()+
  xlab(labs[1]) + ylab(labs[2])+theme(text = element_text(size = 20))
  # geom_text_repel(  
  #   data = to_plot,  
  #   aes(x = PC1, y = PC2, label = rownames, color = age),  
  #   size = 5,  
  #   box.padding = unit(0.35, "lines"),  
  #   point.padding = unit(0.3, "lines")  
  # ) 
p
ggsave("result/Sup_figures/EZHi_RS_PCA.pdf",p,width = 7,height = 6)
