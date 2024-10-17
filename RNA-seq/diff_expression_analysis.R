rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(edgeR)
library(ggplot2)
library(stringr)
tissue <- "FC"
diff_expression_analysis <- function(tissue){
  tab = read.delim(paste0("data/samples/RNA/",tissue,"/combined-chrM.counts"),skip=1)
  rownames(tab) <- tab$Geneid
  tab <- tab[,-1]
  colnames <- colnames(tab)[6:length(tab)]
  pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+).*"
  colnames(tab)[6:length(tab)] <- gsub(pattern, "\\1", colnames(tab)[6:length(tab)] )
  counts <- tab[6:length(tab)]
  group <- read.csv("data/samples/RNA/sample_tissue_info.csv",sep = ',')
  search_table <- read.csv("data/samples/all/RNA_search_table.csv")
  colnames(group)[1] <- "sample_name"
  group <- merge(group,search_table,by="sample_name")
  group <- group[which(group$sample_name %in% colnames(counts)),]
  group$sample_name <- factor(group$sample_name,levels = colnames(counts))
  group <- group[order(group$sample_name),]
  age <- group$Age
  age[which(age=="3m")] <- "young"
  age[which(age=="24m")] <- "old"
  mouse_ID <- group$mouse_ID
  
  colnames(counts) <- paste0(colnames(counts),"-",mouse_ID,"-",age)
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
  write.csv(out,paste0("data/samples/RNA/",tissue,"/diff_expression_gene.csv"))
  # write.csv(out,paste0("data/samples/RNA/DEG_list/",tissue,"_diff_expression_gene_nodup.csv"))
}
tissues <- c("skin","CB","spleen","heart","bladder","tongue","uterus","aorta","thymus","stomach","Hip","FC","BAT","iWAT","muscle","bonemarrow","lung","kidney","liver","testis","colon","cecum","ileum","jejunum","ovary","mammarygland","pancreas")
for(tissue in tissues){
  diff_expression_analysis(tissue)
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
    }else if(tissue_label=="Mammarygland"){
      tissue_label <- "Mammary Gland"
    }else if(tissue_label=="Iwat"){
      tissue_label <- "iWAT"
    }
  }
  return(tissue_label)
}
plot_volcano <- function(tissue){
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene.csv"),row.names = 1)
  colour=setNames(c("blue","grey","red"),c("Down","Stable","Up"))
  p <- ggplot(
    # 数据、映射、颜色
    df, aes(x = logFC, y = -log10(fdr))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = colour) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (fdr)") +
    theme_bw()+
    theme(text = element_text(size = 20))+
    ggtitle(tissue_label_change(tissue))+
    annotate("text", x = min(df$logFC), y = max(-log10(df$fdr)), label = nrow(df[which(df$Significant=="Down"),]), vjust = 5, hjust = 0,colour="blue",size=5)+
    annotate("text", x = max(df$logFC), y = max(-log10(df$fdr)), label = nrow(df[which(df$Significant=="Up"),]), vjust = 5, hjust = 1.5,colour="red",size=5)
  print(p)
  return(p)
}

p_list <- list()

for(i in c(1:length(tissues))){
  p_list[[i]] <- plot_volcano(tissues[i])
}

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_plot <- plot_a_list(p_list,no_of_rows = 4,no_of_cols = 7)
ggsave("result/RNA/all_tissues_diff_expression_gene_volcano.png",combined_plot,width = 40,height = 20,type="cairo")

tissues <- sort(tissues)
diff_gene_number <-data.frame(Var1 = character(),
                              Freq = numeric(),
                              tissue = character(),
                              stringsAsFactors = FALSE)
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene.csv"))
  sig <- data.frame(Var1=c("Up","Stable","Down"),Freq=c(0,0,0))
  t_sig<-as.data.frame(table(df$Significant))
  sig <- merge(sig, t_sig, by="Var1", all.x=TRUE) 
  sig$Freq.x <- ifelse(is.na(sig$Freq.y), sig$Freq.x, sig$Freq.y)  
  colnames(sig)[2] <- "Freq"
  sig <- sig[, -3] 
  sig$tissue <- tissue_label_change(tissue)
  diff_gene_number<-rbind(diff_gene_number,sig)
}

conditions <- c("Down","Up")
for(i in c(1:length(conditions))){
  condition <- conditions[i]
  df <- diff_gene_number[which(diff_gene_number$Var1==condition),]  
  color <- read.table("data/samples/30_distinct_color.txt")
  color <- setNames(color$V1,sort(unique(df$tissue)))
  ggplot(df,mapping = aes(x=Freq,y=tissue,fill = tissue))+
    geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab("RNA")+
    theme(text = element_text(size = 13))+ scale_fill_manual(values = color) +theme(legend.position = "none") + 
    geom_text(aes(label = Freq), position = position_dodge2(width = 0.9), hjust = 0.4, size = 5) + xlim(0,5000)+ggtitle(paste0(condition," gene number"))
}
