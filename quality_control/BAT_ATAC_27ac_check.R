rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)
library(edgeR)
# tab <- read.delim(paste0("data/samples/ATAC/BAT/ATAC/ATAC_1kb_bins.counts"),skip=1)  
# tab <- read.delim(paste0("data/samples/ATAC/BAT/ATAC/ATAC_macs_young_old_narrowpeak.counts"),skip=1)
# search_table <- read.csv("data/samples/all/ATAC_search_table.csv")
tab <- read.delim(paste0("data/samples/BAT/H3K27ac/H3K27ac_1kb_bins.counts"),skip=1)
tab <- read.delim(paste0("data/samples/BAT/H3K27ac/H3K27ac_macs_young_old_narrowpeak.counts"),skip=1)
search_table <- read.csv("data/samples/all/CUTTag_search_table.csv")
pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+).*"
colnames(tab)[7:length(tab)] <- gsub(pattern, "\\1", colnames(tab)[7:length(tab)])
counts <- tab[7:length(tab)]
search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
search_table <- search_table[order(search_table$sample_name),]
mouse_ID <- search_table$mouse_ID
age <- search_table$age
age[which(age=="3m")] <- "young"
age[which(age=="24m")] <- "old"
colnames(counts) <- paste0(colnames(counts),"-",age,"-",mouse_ID)

y= DGEList(counts=counts[,c(1:4)],group = age[1:4])
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y$samples$group <- factor(y$samples$group,c("young","old"))
y <- calcNormFactors(y)
design <- model.matrix(~group, y$samples)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 2)
t_tab<-tab[keep,]
out_previous = cbind(t_tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
            "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
            "LogFC.old-young"=lrt$table$logFC)

y= DGEList(counts=counts[,c(5:8)],group = age[5:8])
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y$samples$group <- factor(y$samples$group,c("young","old"))
y <- calcNormFactors(y)
design <- model.matrix(~group, y$samples)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)
lrt = glmLRT(fit_tag, coef = 2)
t_tab<-tab[keep,]
out_new = cbind(t_tab[,1:6],cpm(y),logCPM=lrt$table$logCPM,bcv=sqrt(fit_tag$dispersion),
                     "PValue.old-young"=lrt$table$PValue,"FDR.old-young"= p.adjust(lrt$table$PValue,method="BH"),
                     "LogFC.old-young"=lrt$table$logFC)

to_plot <- merge(out_previous[,c("Geneid","LogFC.old-young")],out_new[,c("Geneid","LogFC.old-young")],by="Geneid")
colnames(to_plot)[2:3] <- c("previous","new")
ggplot()+
  geom_point(data=to_plot, mapping=aes(previous,new),color = "grey",alpha=0.5) +  
  geom_point(data=to_plot[which(to_plot[,2]>0 & to_plot[,3]>0),], mapping=aes(previous,new),color = "#00b8a9") +
  geom_point(data=to_plot[which(to_plot[,2]<0 & to_plot[,3]<0),], mapping=aes(previous,new),color = "#ff9a00") +
  labs(x="previous",
       y="new") +
  theme_bw()+theme(text = element_text(size = 18))+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggtitle("BAT H3K27ac") +
  annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]>0 & to_plot[,3]>0),])),x = Inf, y = Inf,hjust = 1.1, vjust = 1.2,colour="#00b8a9",size=5)+
  annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]<0 & to_plot[,3]>0),])),x = -Inf, y = Inf,hjust = -0.1, vjust = 1.2,colour="#ff9a00",size=5)+
  annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]<0 & to_plot[,3]<0),])),x = -Inf, y = -Inf,hjust = -0.1, vjust = -1.2,colour="#f6416c",size=5)+
  annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]>0 & to_plot[,3]<0),])),x = Inf, y = -Inf,hjust = 1.1, vjust = -1.2,colour="#48466d",size=5)
par(cex.lab = 1.5, cex.axis = 1.2)
smoothScatter(as.numeric(to_plot[,3]) ~ as.numeric(to_plot[,2]),xlab = colnames(to_plot)[2],ylab = colnames(to_plot)[3],main = "BAT H3K27ac")
abline(a = 0, b = 0, col = "red", lty = 2)
abline(v = 0, col = "red", lty = 2)

# peak <- read.table("data/samples/ATAC/BAT/ATAC/bed/ATAC_1kb_in_young_old_merge_macs_narrowpeak.bed")

to_plot <- to_plot[which(to_plot$Geneid %in% peak$V4),]
ggplot()+
  geom_point(data=to_plot, mapping=aes(previous,new),color = "grey",alpha=0.5) +  
  geom_point(data=to_plot[which(to_plot[,2]>0 & to_plot[,3]>0),], mapping=aes(previous,new),color = "#00b8a9") +
  geom_point(data=to_plot[which(to_plot[,2]<0 & to_plot[,3]<0),], mapping=aes(previous,new),color = "#ff9a00") +
  labs(x="previous",
       y="new") +
  theme_bw()+theme(text = element_text(size = 18))+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggtitle("BAT ATAC") +
  annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]>0 & to_plot[,3]>0),])),x = Inf, y = Inf,hjust = 1.1, vjust = 1.2,colour="#00b8a9",size=5)+
  annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]<0 & to_plot[,3]>0),])),x = -Inf, y = Inf,hjust = -0.1, vjust = 1.2,colour="#ff9a00",size=5)+
  annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]<0 & to_plot[,3]<0),])),x = -Inf, y = -Inf,hjust = -0.1, vjust = -1.2,colour="#f6416c",size=5)+
  annotate("text",label = paste0(nrow(to_plot[which(to_plot[,2]>0 & to_plot[,3]<0),])),x = Inf, y = -Inf,hjust = 1.1, vjust = -1.2,colour="#48466d",size=5)
par(cex.lab = 1.5, cex.axis = 1.2)
smoothScatter(as.numeric(to_plot[,3]) ~ as.numeric(to_plot[,2]),xlab = colnames(to_plot)[2],ylab = colnames(to_plot)[3],main = "BAT ATAC")
abline(a = 0, b = 0, col = "red", lty = 2)
abline(v = 0, col = "red", lty = 2)
