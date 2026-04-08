rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
tissues <- c("FC","Hip","CB")
aging_markers <- c("Apod","Il33","Neat1","Gatm","Ermn","Plp1","Gstp1","Rps21",
                   "App","Cryab","Rpl38","Rpl39","Rpl23a","Pisd","Rpl6","Ptgds",
                   "Rps29","Actb","Tpt1","Malat1","Hsp90aa1","Hsp90ab1","Clu",
                   "Cst3","Mt1","mt-Nd1","Sparcl1","Mt2","Gpr37l1","Cpe","Aldoc",
                   "Ttyh1","Slc1a2","Atp1a2","Mt3")
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
colors <- c("blue","grey","red")
colors <- setNames(colors,c("Down","Stable","Up"))
for(tissue in tissues){
  df <- read.csv(paste0("data/samples/RNA/",tissue,"/diff_expression_gene_nodup.csv"))
  df <- df[which(df$X %in% aging_markers),]
  df<- df[order(df$logFC),]
  df$X <- factor(df$X,levels=df$X)
  ggplot(df,mapping = aes(x=logFC,y=X,fill = Significant))+
    geom_bar(stat = "identity", position = position_dodge2())+theme_bw()+ylab("")+xlab("log2FC")+ggtitle(paste0("Ximerakis aging makers ",tissue_label_change(tissue)))+
    theme(text = element_text(size = 18))+ scale_fill_manual(values = colors) + xlim((min(df$logFC-0.5)),(max(df$logFC+0.5)))+guides(fill= guide_legend(title = ""))
}
