rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(clusterProfiler)
library(stringr)
tissues <- c("skin","CB","spleen","heart","bladder","tongue","uterus","aorta","thymus","stomach","Hip","FC","BAT","iWAT","muscle","bonemarrow","lung","kidney","liver","testis","colon","cecum","ileum","jejunum")
# tissues <- c("bladder","tongue","uterus","aorta","thymus","stomach","Hip","FC","BAT","iWAT","muscle","bonemarrow","lung","kidney","liver","testis","colon","cecum","ileum","jejunum")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
GO_database <- 'org.Mm.eg.db'
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
    }
  }
  return(tissue_label)
} 
dir.create("result/RNA/GO/plot")
dir.create("result/RNA/GO/table")
GO_pathway_analysis <- function(tissue){
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_nodup.csv"),row.names = 1)
  genelist_up <- bitr(rownames(df)[which(df$Significant=="Up")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_up_GO <- enrichGO( genelist_up$ENTREZID,#GO富集分析
                              OrgDb = GO_database,
                              keyType = "ENTREZID",#设定读取的gene ID类型
                              ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                              pvalueCutoff = 0.05,#设定p值阈值
                              qvalueCutoff = 0.05,#设定q值阈值
                              readable = T)
  if(nrow(genelist_up_GO)>0){
    p <- barplot(genelist_up_GO,title = paste0(tissue_label_change(tissue)," Increased gene GO pathway"),label_format = 50)
    ggsave(paste0("result/RNA/GO/plot/",tissue_label_change(tissue),"_increased_gene_GO.png"),p,width = 10,height = 5,type="cairo")
  }
  GO_table <- genelist_up_GO@result
  write.csv(GO_table,paste0("result/RNA/GO/table/",tissue_label_change(tissue),"_increased_gene_GO.csv"))

  genelist_down <- bitr(rownames(df)[which(df$Significant=="Down")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genelist_down_GO <- enrichGO( genelist_down$ENTREZID,#GO富集分析
                                OrgDb = GO_database,
                                keyType = "ENTREZID",#设定读取的gene ID类型
                                ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                                pvalueCutoff = 0.05,#设定p值阈值
                                qvalueCutoff = 0.05,#设定q值阈值
                                readable = T)
  if(nrow(genelist_down_GO)>0){
    p <- barplot(genelist_down_GO,title = paste0(tissue_label_change(tissue)," Decreased gene GO pathway"),label_format = 50)
    ggsave(paste0("result/RNA/GO/plot/",tissue_label_change(tissue),"_decreased_gene_GO.png"),p,width = 10,height = 5,type="cairo")
  }
  GO_table <- genelist_down_GO@result
  write.csv(GO_table,paste0("result/RNA/GO/table/",tissue_label_change(tissue),"_decreased_gene_GO.csv"))

}

for (tissue in tissues){
  GO_pathway_analysis(tissue)
}

geneup_list <- list()
genedown_list <- list()
for (i in c(1:length(tissues))){
  tissue <- tissues[i]
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_nodup.csv"),row.names = 1)
  gene_up <- bitr(rownames(df)[which(df$Significant=="Up")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  if(length(gene_up)>0){
    geneup_list[[i]] <- gene_up$ENTREZID
    names(geneup_list)[i] <- tissue_label_change(tissue)
  }
  gene_down <- bitr(rownames(df)[which(df$Significant=="Down")],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  if(length(gene_up)>0){
    genedown_list[[i]] <- gene_down$ENTREZID
    names(genedown_list)[i] <- tissue_label_change(tissue)
  }
}
ck <- compareCluster(geneCluster = geneup_list, fun = enrichGO,OrgDb = GO_database, keyType = "ENTREZID",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
p<- dotplot(ck,show=5,label_format = 50)+    
  theme(  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),                           
    axis.title = element_text(size = 14),                           
    plot.title = element_text(size = 16),                          
    legend.text = element_text(size = 12),                         
    legend.title = element_text(size = 14)                        
  ) +labs(x=NULL)
ggsave(paste0("result/RNA/GO/plot/all_tissues_increase_GO_pathway.png"),p,width = 20,height = 20,type="cairo")

ck <- compareCluster(geneCluster = genedown_list, fun = enrichGO,OrgDb = GO_database, keyType = "ENTREZID",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
p<- dotplot(ck,show=5,label_format = 100)+    
  theme(  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),                           
    axis.title = element_text(size = 14),                           
    plot.title = element_text(size = 16),                          
    legend.text = element_text(size = 12),                         
    legend.title = element_text(size = 14)                        
  ) +labs(x=NULL)
ggsave(paste0("result/RNA/GO/plot/all_tissues_decrease_GO_pathway.png"),p,width = 20,height = 20,type="cairo")
