rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library("AnnotationDbi")
library(org.Mm.eg.db)
library(edgeR)
library(ggplot2)
library(ChIPseeker)
library(EnsDb.Mmusculus.v79)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggrepel)
library(limma)
library(ViSEAGO)
library(data.table)  
library(DOSE)
library(enrichplot)
options(bitmapType='cairo')
bin_size <- function(antibody){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    return ("10kb")
  }else{
    return("1kb")
  }  
}
# pattern.extract <- function(query, m) {  
#   Data = lapply(seq_along(query), function(i) {  
#     if (length(stats::na.omit(m[[i]][1])) > 0) {  
#       capture = attr(m[[i]], "capture.start")  
#       capture = substring(query[i], capture, capture +   
#                             attr(m[[i]], "capture.length") - 1)  
#       data.table(t(capture))  
#     } else {  
#       NA  
#     }  
#   })  
#   Data <- rbindlist(Data)  
#   colnames(Data) <- attr(m[[1]], "capture.name")  
#   Data[c("", "\t"), `:=`("CommonName", "NA"), on = "CommonName"]  
#   return(Data)  
# }  
# 
# temp <- ("/storage/zhangyanxiaoLab/suzhuojie/ref_data/EntrezGene/gene2go.gz")
# system(paste("gunzip",temp)) 
# gene2go = fread(sub("\\.gz", "", temp), verbose = FALSE, 
#                 showProgress = FALSE)
# gene2go <- unique(gene2go[, c(seq_len(4), 8), with = FALSE])
# colnames(gene2go) <- c("taxid", "gene_id", "GOID", "evidence", 
#                        "category")
# gene2go[, `:=`(taxid = as.character(gene2go$taxid), gene_id = as.character(gene2go$gene_id))]
# core = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?version=2.0&db=taxonomy"  
# taxid  <- unique(gene2go$taxid)
# for(i in c(1:length(taxid))){
#   query <- paste(core, "&id=", paste(taxid[i], collapse = ","), sep = "")
#   query = paste(scan(query, what = "", sep = "\n", quiet = TRUE), collapse = "")  
#   query <- substring(query, unlist(gregexpr("<DocumentSummary ", query)), unlist(gregexpr("</DocumentSummary>", query)))  
#   pattern = c("<DocumentSummary uid=\"(?<taxid>[[:digit:]]*)", ".*<ScientificName>(?<ScientificName>.*)</ScientificName>", "\t<CommonName>(?<CommonName>.*)</CommonName>.*")
#   m = gregexpr(paste(pattern, collapse = ""), query, perl = TRUE) 
#   result <- pattern.extract(query, m) 
#   if(i == 1) {
#     taxon <- result
#   }else{
#     taxon <- rbind(taxon, result)
#   }
# }
# EntrezGene <- new("genomic_ressource", db = "EntrezGene", stamp = as.character(Sys.time()), 
#     data = gene2go, organisms = taxon)
# saveRDS(EntrezGene,"~/ref_data/EntrezGene/EntrezGene.rds")
# ViSEAGO::available_organisms(EntrezGene)
# myGENE2GO<-ViSEAGO::annotate(
#   "10090",
#   EntrezGene
# )

antibody<-"H3K27ac"
tissues = c("brain","liver","testis","colon","kidney","lung","spleen","muscle","Hip","cecum","bonemarrow","ileum","heart","thymus")
GO_database <- 'org.Mm.eg.db'
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
geneup_list <- list()
genedown_list <- list()
for(i in c(1:length(tissues))){
  tissue <- tissues[i]
  out <- read.csv(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size(antibody),"_bins_diff_after_remove_batch_effect.csv"))
  peak <- GRanges(seqnames = out$Chr,   
                  ranges = IRanges(start = out$Start, end = out$End))
  peak_anno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  peak_anno <- unique(as.data.frame(peak_anno))
  genelist <- bitr(peak_anno$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  background <- genelist$ENTREZID
  out_up <- out[which(out$Significant_bar=="Up"),]
  up_peak <- GRanges(seqnames = out_up$Chr,   
                     ranges = IRanges(start = out_up$Start, end = out_up$End))
  up_peak_anno <- annotatePeak(up_peak, tssRegion=c(-3000, 3000),
                               TxDb=txdb, annoDb="org.Mm.eg.db")
  up_peak_anno <- unique(as.data.frame(up_peak_anno))
  genelist_up <- bitr(up_peak_anno$SYMBOL[which(str_detect(up_peak_anno$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  geneup_list[[i]] <- genelist_up$ENTREZID
  names(geneup_list)[i] <- tissue
  out_down <- out[which(out$Significant_bar=="Down"),]
  down_peak <- GRanges(seqnames = out_down$Chr,   
                     ranges = IRanges(start = out_down$Start, end = out_down$End))
  down_peak_anno <- annotatePeak(down_peak, tssRegion=c(-3000, 3000),
                               TxDb=txdb, annoDb="org.Mm.eg.db")
  down_peak_anno <- unique(as.data.frame(down_peak_anno))
  genelist_down <- bitr(down_peak_anno$SYMBOL[which(str_detect(down_peak_anno$annotation,"Promoter"))],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  genedown_list[[i]] <- genelist_down$ENTREZID
  names(genedown_list)[i] <- tissue

}
ck <- compareCluster(geneCluster = geneup_list, fun = enrichGO,OrgDb = GO_database, keyType = "ENTREZID",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
p<- dotplot(ck,show=3,label_format = 100)+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +labs(x=NULL)
ggsave(paste0("result/all/diff/",antibody,"/1kb_increase_promoter_Go_dotplot.png"),p,width = 13,height = 10)
ck <- compareCluster(geneCluster = genedown_list, fun = enrichGO,OrgDb = GO_database, keyType = "ENTREZID",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
p<- dotplot(ck,show=3,label_format = 100)+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +labs(x=NULL)
ggsave(paste0("result/all/diff/",antibody,"/1kb_decrease_promoter_Go_dotplot.png"),p,width = 13,height = 10)

# cnetplot(ck)
# Input <- list()
# for(i in c(1:length(tissues))){
#   Input[[i]] <- c(tissues[i],paste0(tissues[i],"_classic"))
#   names(Input)[i] <- tissues[i]
# }
# BP_sResults<-ViSEAGO::merge_enrich_terms(
#   Input=Input
# )
# myGOs<-ViSEAGO::build_GO_SS(
#   gene2GO=myGENE2GO,
#   enrich_GO_terms=BP_sResults
# )
# myGOs<-ViSEAGO::compute_SS_distances(
#   myGOs,
#   distance="Wang"
# )
# ViSEAGO::MDSplot(myGOs)
# Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#   myGOs,
#   showIC=TRUE,
#   showGOlabels=TRUE,
#   GO.tree=list(
#     tree=list(
#       distance="Wang",
#       aggreg.method="ward.D2"
#     ),
#     cut=list(
#       dynamic=list(
#         pamStage=TRUE,
#         pamRespectsDendro=TRUE,
#         deepSplit=2,
#         minClusterSize =2
#       )
#     )
#   ),
#   samples.tree=NULL
# )
# ViSEAGO::show_heatmap(
#   Wang_clusters_wardD2,
#   "GOterms"
# )
# ViSEAGO::show_table(
#   Wang_clusters_wardD2,
#   "cluster_heatmap_Wang_wardD2.xls"
# )
# Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(
#     Wang_clusters_wardD2,
#     distance=c("max", "avg","rcmax", "BMA")
# )
# ViSEAGO::MDSplot(
#   Wang_clusters_wardD2,
#   "GOclusters"
# ) 
# Wang_clusters_wardD2<-ViSEAGO::GOclusters_heatmap(
#   Wang_clusters_wardD2,
#   tree=list(
#     distance="BMA",
#     aggreg.method="ward.D2"
#   )
# )
# ViSEAGO::show_heatmap(
#   Wang_clusters_wardD2,
#   "GOclusters"
# )
# 
