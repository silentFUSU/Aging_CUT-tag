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
tissue <- "skin"
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

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols,guides = "collect")
}
tissue_label_change <- function(tissue){
  if(tissue=="FC"){
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
    }else if (tissue_label == "Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
} 

TE_diff_analysis <- function(tissue){
  tab <- read.delim(paste0("data/samples/RNA/",tissue,"/TEcount/combined.cntTable"),row.names = 1)
  counts <- tab
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SRR[0-9]+|HM[0-9]+).*"
  colnames(counts) <- gsub(pattern, "\\1", colnames(counts))
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
  search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
  search_table <- search_table[order(search_table$sample_name),]
  age <- search_table$age
  mouse_ID <- search_table$mouse_ID
  colnames(counts) <- paste0(colnames(counts),"-",mouse_ID,"-",age)
  age[which(age=="3m")] <- "young"
  age[which(age=="24m")] <- "old"
  
  y= DGEList(counts=counts,group=age)
  # keep = which(rowSums(cpm(y)>1)>=2)
  keep = which(rowSums(cpm(y)>0)>=2)
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
    head(5)  
  split_names <- strsplit(rownames(top_TE_Down), ":")  
  split_df <- do.call(rbind, split_names) 
  top_TE_Down <- cbind(top_TE_Down, split_df)
  colnames(top_TE_Down)[ncol(top_TE_Down)-2] <- "gene"
  
  top_TE_Up <- out %>%  
    filter(Significant == "TE-Up") %>%  
    arrange(fdr) %>%  
    head(5)  
  split_names <- strsplit(rownames(top_TE_Up), ":")  
  split_df <- do.call(rbind, split_names) 
  top_TE_Up <- cbind(top_TE_Up, split_df)
  colnames(top_TE_Up)[ncol(top_TE_Up)-2] <- "gene"
  out <- out[which(!grepl("^ENSMUSG", rownames(out))), ]
  # write.csv(out,paste0("data/samples/RNA/",tissue,"/diff_expression_TE.csv"))
  write.csv(out,paste0("data/samples/RNA/",tissue,"/diff_expression_TE_change_filter_bar.csv"))
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
    ggtitle(tissue_label_change(tissue))+
    annotate("text", x = min(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="TE-Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(out$logFC), y = max(-log10(out$fdr)), label = nrow(out[which(out$Significant=="TE-Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  p <- p+geom_text_repel(data=top_TE_Down, mapping=aes(x=logFC, y=-log10(fdr), label=gene), vjust=-1, size=4) +  
    geom_text_repel(data=top_TE_Up, mapping=aes(x=logFC, y=-log10(fdr), label=gene), vjust=-1, size=4)  
  return(p)
  
} 

tissues <- sort(c("skin","CB","spleen","heart","bladder","tongue","uterus","aorta","thymus","stomach","Hip","brain","BAT","iWAT","muscle","bonemarrow","lung","kidney","liver","testis","colon","cecum","ileum","jejunum","pancreas","ovary","mammarygland"))
p_list <- list()
for (i in c(1:length(tissues))){
  p_list[[i]] <- TE_diff_analysis(tissues[i])
}
combined_plot <- plot_a_list(p_list, no_of_rows = 4,no_of_cols = 7)
ggsave("result/RNA/TE/TE_volcano_all_tissues_change_filter_bar.png",combined_plot,width = 30,height = 20,type="cairo")


diff_TE_number <-data.frame(Var1 = character(),
                              Freq = numeric(),
                              tissue = character(),
                              stringsAsFactors = FALSE)
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_TE.csv"))
  sig <- data.frame(Var1=c("TE-Up","Stable","TE-Down"),Freq=c(0,0,0))
  t_sig<-as.data.frame(table(df$Significant))
  sig <- merge(sig, t_sig, by="Var1", all.x=TRUE) 
  sig$Freq.x <- ifelse(is.na(sig$Freq.y), sig$Freq.x, sig$Freq.y)  
  colnames(sig)[2] <- "Freq"
  sig <- sig[, -3] 
  sig$tissue <- tissue_label_change(tissue)
  diff_TE_number<-rbind(diff_TE_number,sig)
}

conditions <- c("TE-Down","TE-Up")
for(i in c(1:length(conditions))){
  condition <- conditions[i]
  df <- diff_TE_number[which(diff_TE_number$Var1==condition),]  
  color <- read.table("data/samples/30_distinct_color.txt")
  color <- setNames(color$V1,sort(unique(df$tissue)))
  p <- ggplot(df,mapping = aes(x=Freq,y=tissue,fill = tissue))+
    geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab("TE")+
    theme(text = element_text(size = 13))+ scale_fill_manual(values = color) +theme(legend.position = "none") + 
    geom_text(aes(label = Freq), position = position_dodge2(width = 0.9), hjust = 0.4, size = 5) + xlim(0,700)+ggtitle(paste0(condition," gene number"))
  print(p)
  }
