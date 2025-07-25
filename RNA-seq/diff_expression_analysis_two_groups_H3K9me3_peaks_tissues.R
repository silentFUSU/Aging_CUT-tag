rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
library(ggrepel)
tab = read.delim(paste0("data/samples/RNA/all_tissues_combined-chrM.counts"),skip=1)

rownames(tab) <- tab$Geneid
tab <- tab[,-1]
colnames <- colnames(tab)[6:length(tab)]
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|HM[0-9]+).*"
colnames(tab)[6:length(tab)] <- gsub(pattern, "\\1", colnames(tab)[6:length(tab)] )
counts <- tab[6:length(tab)]

group <- read.csv("data/samples/RNA/sample_tissue_info.csv",sep = ',')
search_table <- read.csv("data/samples/all/RNA_search_table.csv")
colnames(group)[1] <- "sample_name"
group <- merge(group,search_table,by="sample_name")

tissues1 <- c("Cortex","Hippocampus","Skin","Heart","Aorta","Muscle","Cerebellum","BAT","Lung")
tissues2 <- c("Ileum","Spleen","Thymus","Liver","Cecum","Colon","Jejunum","Bone Marrow","Stomach","Kidney","Pancreas","iWAT","Bladder","Testis","Tongue")
group <- group[which(group$tissue %in% c(tissues1,tissues2)),]

counts <- counts[,group$sample_name]
group <- group[which(group$sample_name %in% colnames(counts)),]
group$sample_name <- factor(group$sample_name,levels = colnames(counts))
group <- group[order(group$sample_name),]
group$condition <- "group1"
group$condition[which(group$tissue%in%tissues2)] <- "group2"

condition <- group$condition
y= DGEList(counts=counts,group=condition)
keep = which(rowSums(cpm(y)>0)>=2)
y = y[keep,]
y$samples$group <- factor(y$samples$group, levels=c("group1","group2"))
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

top_genes <- out[which(out$Significant=="Up"),] %>%  
  arrange(desc(logFC)) %>%  
  head(10)  
bottom_genes <- out[which(out$Significant=="Down"),] %>%  
  arrange(logFC) %>%  
  head(10)  
highlight_genes <- rbind(top_genes, bottom_genes) 
highlight_genes$Gene <- rownames(highlight_genes)
colour=setNames(c("blue","grey","red"),c("Down","Stable","Up"))
ggplot(
  out, aes(x = logFC, y = -log10(fdr))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = colour) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (fdr)") +
  theme_bw()+
  theme(text = element_text(size = 20))+
  ggtitle("group2 vs group1")+
  annotate("text", x = min(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
  annotate("text", x = max(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)+
  geom_text_repel(data = highlight_genes, aes(logFC,-log10(fdr), label = Gene), max.overlaps=100,
                  size = 5, 
                  nudge_y = 0.2)

cpm <- as.data.frame(cpm(y))
cpm <- cpm[which(rownames(cpm) %in% rownames(out)[which(out$Significant!="Stable")]),]
annotation <- group[,c("sample_name","condition","tissue")]
rownames(annotation) <- annotation$sample_name
annotation <- annotation[,-1,drop=F]

breaks <- c(seq(-4, -1.1, length.out = 40), seq(-1, 1, length.out = 20), seq(1.1, 4, length.out = 40))
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap::pheatmap(cpm,show_rownames = F,scale="row",annotation_col = annotation,breaks = breaks,color = color_palette)
