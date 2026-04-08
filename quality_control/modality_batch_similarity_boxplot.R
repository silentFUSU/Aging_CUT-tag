rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(ggplot2)

tissue_label_change <- function(tissue){
  if(tissue=="brain"){
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
    }
  }
  return(tissue_label)
}

### CUT&Tag
antibodys <- c("H3K27me3","H3K9me3","H3K36me3","H3K27ac","H3K4me1","H3K4me3")
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder","CB","jejunum","uterus","ovary")
p_list <- list()
for(antibody in antibodys){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  antibody_correlation <- data.frame()
  for(tissue in tissues){
    search_table <- read.csv("data/samples/all/CUTTag_search_table_used_in_diff_batch.csv")
    tab <- read.table(paste0("data/samples/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),header = T)
    if(tissue %in% c("mammarygland","uterus","ovary","MEF")){
      tab <- tab[which(!tab$Chr %in% c("chrY")),]
    }
    rownames(tab) <- tab$Geneid
    counts = tab[,c(7:ncol(tab))]
    rownames(counts)= tab$Geneid

    pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|XM[0-9]+|DYQ[0-9]+).*"
    colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
    search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
    counts <- counts[,search_table$sample_name]
    search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
    search_table$age <- factor(search_table$age, levels = c("3m","24m"))
    search_table <- search_table[order(search_table$age),]
    search_table$batch <- rep(paste0("batch",c(1:(nrow(search_table)/2))),2)
    search_table <- search_table[which(search_table$batch %in% c("batch1","batch2")),]
    counts <- counts[,search_table$sample_name]
    y= DGEList(counts=counts)
    keep = which(rowSums(cpm(y)>1)>=2)
    y = y[keep,]
    cpm <- as.data.frame(cpm(y, log = TRUE))
    vars <- apply(cpm, 1, var, na.rm = TRUE)
    cpm <- cpm[order(vars, decreasing = TRUE)[1:100000], ]
    cpm <- as.data.frame(limma::removeBatchEffect(cpm, batch = search_table$batch))
    search_table$condition <- paste0(search_table$age,"_",search_table$batch)
    search_table$condition[which(search_table$condition=="3m_batch1")] <- "Y1"
    search_table$condition[which(search_table$condition=="3m_batch2")] <- "Y2"
    search_table$condition[which(search_table$condition=="24m_batch1")] <- "O1"
    search_table$condition[which(search_table$condition=="24m_batch2")] <- "O2"
    pair_mat <- combn(search_table$condition, 2)
    pair_df  <- as.data.frame(t(pair_mat))
    
    
    tissue_correaltion <- data.frame()
    for(i in 1:nrow(pair_df)){
      condition1 <- pair_df$V1[i]
      condition2 <- pair_df$V2[i]
      t_cpm <- cpm[,as.character(c(search_table$sample_name[which(search_table$condition==condition1)],search_table$sample_name[which(search_table$condition==condition2)]))]
      cor <- cor.test(t_cpm[,1],t_cpm[,2])
      t_tissue_correlation <- data.frame(tissue=tissue_label_change(tissue),cor=as.numeric(cor$estimate[1]),pair=paste0(condition1,"-",condition2))
      tissue_correaltion <- rbind(tissue_correaltion,t_tissue_correlation)
    }
    antibody_correlation <- rbind(antibody_correlation,tissue_correaltion)
  }
  antibody_correlation$pair <- factor(antibody_correlation$pair,levels=c("Y1-Y2","O1-O2","Y1-O1","Y2-O2","Y1-O2","Y2-O1"))
  antibody_correlation$pair2 <- "Y-O"
  antibody_correlation$pair2[which(antibody_correlation$pair=="Y1-Y2")] <- "Y1-Y2"
  antibody_correlation$pair2[which(antibody_correlation$pair=="O1-O2")] <-"O1-O2"
  antibody_correlation$pair2 <- factor(antibody_correlation$pair2,levels=c("Y1-Y2","O1-O2","Y-O"))
  # color <- setNames(c("#f39b7f","#f39b7f","#4dbbd5","#4dbbd5","#4dbbd5","#4dbbd5"),c("Y1-Y2","O1-O2","Y1-O1","Y1-O2","Y2-O1","Y2-O2"))
  color <- setNames(c("#f39b7f","#4dbbd5","grey"),c("Y1-Y2","O1-O2","Y-O")) 
  p_list[[antibody]] <- ggplot(antibody_correlation, aes(x = pair2, y = cor,fill=pair2)) +
    geom_boxplot(outliers =F) +
    scale_fill_manual(values=color)+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab(NULL) + ylab("Correlation")+
    ggtitle(antibody)
}
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}
combined_plot <- plot_a_list(p_list,2,3)
for(antibody in antibodys){
  p <- p_list[[antibody]]
  ggsave(paste0("result/Sup_figures/",antibody,"_correlation_boxplot.pdf"),p,width = 6,height = 8)
  }

### ATAC
antibodys <- c("ATAC")
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder","CB","jejunum","uterus","ovary")
for(antibody in antibodys){
  if(antibody %in% c("H3K27me3","H3K9me3","H3K36me3")){
    bin_size <- "10kb"
  }else{
    bin_size <- "1kb"
  }
  antibody_correlation <- data.frame()
  for(tissue in tissues){
    search_table <- read.csv("data/samples/all/ATAC_search_table_batch.csv")
    tab <- read.table(paste0("data/samples/ATAC/",tissue,"/",antibody,"/",antibody,"_",bin_size,"_bins.counts"),header = T)
    if(tissue %in% c("mammarygland","uterus","ovary","MEF")){
      tab <- tab[which(!tab$Chr %in% c("chrY")),]
    }
    rownames(tab) <- tab$Geneid
    counts = tab[,c(7:ncol(tab))]
    rownames(counts)= tab$Geneid
    
    pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|XM[0-9]+|DYQ[0-9]+).*"
    colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
    search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
    counts <- counts[,search_table$sample_name]
    search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
    search_table$age <- factor(search_table$age, levels = c("3m","24m"))
    search_table <- search_table[order(search_table$age),]
    search_table$batch <- rep(paste0("batch",c(1:(nrow(search_table)/2))),2)
    search_table <- search_table[which(search_table$batch %in% c("batch1","batch2")),]
    counts <- counts[,search_table$sample_name]
    y= DGEList(counts=counts)
    keep = which(rowSums(cpm(y)>1)>=2)
    y = y[keep,]
    cpm <- as.data.frame(cpm(y, log = TRUE))
    vars <- apply(cpm, 1, var, na.rm = TRUE)
    cpm <- cpm[order(vars, decreasing = TRUE)[1:100000], ]
    cpm <- as.data.frame(limma::removeBatchEffect(cpm, batch = search_table$batch))
    
    
    search_table$condition <- paste0(search_table$age,"_",search_table$batch)
    search_table$condition[which(search_table$condition=="3m_batch1")] <- "Y1"
    search_table$condition[which(search_table$condition=="3m_batch2")] <- "Y2"
    search_table$condition[which(search_table$condition=="24m_batch1")] <- "O1"
    search_table$condition[which(search_table$condition=="24m_batch2")] <- "O2"
    pair_mat <- combn(search_table$condition, 2)
    pair_df  <- as.data.frame(t(pair_mat))
    
    
    tissue_correaltion <- data.frame()
    for(i in 1:nrow(pair_df)){
      condition1 <- pair_df$V1[i]
      condition2 <- pair_df$V2[i]
      t_cpm <- cpm[,as.character(c(search_table$sample_name[which(search_table$condition==condition1)],search_table$sample_name[which(search_table$condition==condition2)]))]
      cor <- cor.test(t_cpm[,1],t_cpm[,2])
      t_tissue_correlation <- data.frame(tissue=tissue_label_change(tissue),cor=as.numeric(cor$estimate[1]),pair=paste0(condition1,"-",condition2))
      tissue_correaltion <- rbind(tissue_correaltion,t_tissue_correlation)
    }
    antibody_correlation <- rbind(antibody_correlation,tissue_correaltion)
  }
  antibody_correlation$pair <- factor(antibody_correlation$pair,levels=c("Y1-Y2","O1-O2","Y1-O1","Y2-O2","Y1-O2","Y2-O1"))
  antibody_correlation$pair2 <- "Y-O"
  antibody_correlation$pair2[which(antibody_correlation$pair=="Y1-Y2")] <- "Y1-Y2"
  antibody_correlation$pair2[which(antibody_correlation$pair=="O1-O2")] <-"O1-O2"
  antibody_correlation$pair2 <- factor(antibody_correlation$pair2,levels=c("Y1-Y2","O1-O2","Y-O"))
  # color <- setNames(c("#f39b7f","#f39b7f","#4dbbd5","#4dbbd5","#4dbbd5","#4dbbd5"),c("Y1-Y2","O1-O2","Y1-O1","Y1-O2","Y2-O1","Y2-O2"))
  color <- setNames(c("#f39b7f","#4dbbd5","grey"),c("Y1-Y2","O1-O2","Y-O")) 
  p_list[[antibody]] <- ggplot(antibody_correlation, aes(x = pair2, y = cor,fill=pair2)) +
    geom_boxplot(outliers =F) +
    scale_fill_manual(values=color)+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab(NULL) + ylab("Correlation")+
    ggtitle(antibody)
}
for(antibody in antibodys){
  p <- p_list[[antibody]]
  ggsave(paste0("result/Sup_figures/",antibody,"_correlation_boxplot.pdf"),p,width = 6,height = 8)
}

#### RNA
antibodys <- c("RNA")
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder","CB","jejunum","uterus","ovary")
for(antibody in antibodys){
  antibody_correlation <- data.frame()
  for(tissue in tissues){
    search_table <- read.csv("data/samples/all/RNA_search_table.csv")
    tab <- read.table(paste0("data/samples/RNA/",tissue,"/combined-chrM.counts"),header = T)
    if(tissue %in% c("mammarygland","uterus","ovary","MEF")){
      tab <- tab[-grep("chrY", tab$Chr), ]
    }
    rownames(tab) <- tab$Geneid
    counts = tab[,c(7:ncol(tab))]
    rownames(counts)= tab$Geneid
    
    pattern <- ".*bam\\.(LLX[0-9]+|CKJ[0-9]+|SZJ[0-9]+|HJC[0-9]+|HJC_[0-9]+|NTY[0-9]+|XM[0-9]+|DYQ[0-9]+|HM[0-9]+).*"
    colnames(counts) <-  gsub(pattern, "\\1",colnames(counts))
    search_table <- search_table[which(search_table$sample_name %in% colnames(counts)),]
    counts <- counts[,search_table$sample_name]
    search_table$sample_name <- factor(search_table$sample_name, levels = colnames(counts))
    search_table$age <- factor(search_table$age, levels = c("3m","24m"))
    search_table <- search_table[order(search_table$age),]
    if(nrow(search_table==5)){
      search_table<-search_table[1:4,]
    }
    search_table$batch <- rep(paste0("batch",c(1:(nrow(search_table)/2))),2)
    search_table <- search_table[which(search_table$batch %in% c("batch1","batch2")),]
    counts <- counts[,search_table$sample_name]
    y= DGEList(counts=counts)
    keep = which(rowSums(cpm(y)>1)>=2)
    y = y[keep,]
    cpm <- as.data.frame(cpm(y, log = TRUE))
    vars <- apply(cpm, 1, var, na.rm = TRUE)
    cpm <- cpm[order(vars, decreasing = TRUE)[1:10000], ]
    cpm <- as.data.frame(limma::removeBatchEffect(cpm, batch = search_table$batch))
    
    
    search_table$condition <- paste0(search_table$age,"_",search_table$batch)
    search_table$condition[which(search_table$condition=="3m_batch1")] <- "Y1"
    search_table$condition[which(search_table$condition=="3m_batch2")] <- "Y2"
    search_table$condition[which(search_table$condition=="24m_batch1")] <- "O1"
    search_table$condition[which(search_table$condition=="24m_batch2")] <- "O2"
    pair_mat <- combn(search_table$condition, 2)
    pair_df  <- as.data.frame(t(pair_mat))
    
    
    tissue_correaltion <- data.frame()
    for(i in 1:nrow(pair_df)){
      condition1 <- pair_df$V1[i]
      condition2 <- pair_df$V2[i]
      t_cpm <- cpm[,as.character(c(search_table$sample_name[which(search_table$condition==condition1)],search_table$sample_name[which(search_table$condition==condition2)]))]
      cor <- cor.test(t_cpm[,1],t_cpm[,2])
      t_tissue_correlation <- data.frame(tissue=tissue_label_change(tissue),cor=as.numeric(cor$estimate[1]),pair=paste0(condition1,"-",condition2))
      tissue_correaltion <- rbind(tissue_correaltion,t_tissue_correlation)
    }
    antibody_correlation <- rbind(antibody_correlation,tissue_correaltion)
  }
  antibody_correlation$pair <- factor(antibody_correlation$pair,levels=c("Y1-Y2","O1-O2","Y1-O1","Y2-O2","Y1-O2","Y2-O1"))
  antibody_correlation$pair2 <- "Y-O"
  antibody_correlation$pair2[which(antibody_correlation$pair=="Y1-Y2")] <- "Y1-Y2"
  antibody_correlation$pair2[which(antibody_correlation$pair=="O1-O2")] <-"O1-O2"
  antibody_correlation$pair2 <- factor(antibody_correlation$pair2,levels=c("Y1-Y2","O1-O2","Y-O"))
  # color <- setNames(c("#f39b7f","#f39b7f","#4dbbd5","#4dbbd5","#4dbbd5","#4dbbd5"),c("Y1-Y2","O1-O2","Y1-O1","Y1-O2","Y2-O1","Y2-O2"))
  color <- setNames(c("#f39b7f","#4dbbd5","grey"),c("Y1-Y2","O1-O2","Y-O")) 
  p_list[[antibody]] <- ggplot(antibody_correlation, aes(x = pair2, y = cor,fill=pair2)) +
    geom_boxplot(outliers =F) +
    scale_fill_manual(values=color)+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab(NULL) + ylab("Correlation")+
    ggtitle(antibody)
}
for(antibody in antibodys){
  p <- p_list[[antibody]]
  ggsave(paste0("result/Sup_figures/",antibody,"_correlation_boxplot.pdf"),p,width = 6,height = 8)
}

### HiC
antibodys <- c("HiC")
tissues <- c("brain","CB","kidney","liver","lung","bonemarrow","colon","heart","Hip","mammarygland","stomach","thymus","skin","muscle","cecum","ileum","pancreas","spleen")
resolution <- "50000"
for(antibody in antibodys){
  antibody_correlation <- data.frame()
  for(tissue in tissues){
    search_table <- read.csv("data/samples/all/HiC_search_table.csv")
    search_table <- search_table[which(search_table$tissue==tissue),]
    df_pca <- data.frame()
    for(i in c(1:length(search_table$sample_name))){
      sample <- search_table$sample_name[i]
      df <- read.table(paste0("data/samples/HiC/",tissue,"/compartment/homer_compartment/PC1/",sample,"_",resolution,".PC1.txt"))
      df <- df[which(df$V2 %in% c(paste0("chr",c(1:19,"X")))),]
      df <- df[,c(1,6)]
      colnames(df)[2] <- sample
      if(nrow(df_pca)==0){
        df_pca <- df
      }else{
        df_pca <- merge(df_pca,df,by="V1")
      }
    }
    rownames(df_pca) <- df_pca$V1
    df_pca <- df_pca[,-1]
    vars <- apply(df_pca, 1, var, na.rm = TRUE)
    df_pca <- df_pca[order(vars, decreasing = TRUE)[1:10000], ]
    
    search_table$age <- factor(search_table$age, levels = c("3M","24M"))
    search_table <- search_table[order(search_table$age),]
    search_table$batch <- rep(paste0("batch",c(1:(nrow(search_table)/2))),2)
    search_table <- search_table[which(search_table$batch %in% c("batch1","batch2")),]
    
    search_table$condition <- paste0(search_table$age,"_",search_table$batch)
    search_table$condition[which(search_table$condition=="3M_batch1")] <- "Y1"
    search_table$condition[which(search_table$condition=="3M_batch2")] <- "Y2"
    search_table$condition[which(search_table$condition=="24M_batch1")] <- "O1"
    search_table$condition[which(search_table$condition=="24M_batch2")] <- "O2"
    pair_mat <- combn(search_table$condition, 2)
    pair_df  <- as.data.frame(t(pair_mat))
    
    
    tissue_correaltion <- data.frame()
    for(i in 1:nrow(pair_df)){
      condition1 <- pair_df$V1[i]
      condition2 <- pair_df$V2[i]
      t_df_pca <- df_pca[,as.character(c(search_table$sample_name[which(search_table$condition==condition1)],search_table$sample_name[which(search_table$condition==condition2)]))]
      cor <- cor.test(t_df_pca[,1],t_df_pca[,2])
      t_tissue_correlation <- data.frame(tissue=tissue_label_change(tissue),cor=as.numeric(cor$estimate[1]),pair=paste0(condition1,"-",condition2))
      tissue_correaltion <- rbind(tissue_correaltion,t_tissue_correlation)
    }
    antibody_correlation <- rbind(antibody_correlation,tissue_correaltion)
  }
  antibody_correlation$pair <- factor(antibody_correlation$pair,levels=c("Y1-Y2","O1-O2","Y1-O1","Y2-O2","Y1-O2","Y2-O1"))
  antibody_correlation$pair2 <- "Y-O"
  antibody_correlation$pair2[which(antibody_correlation$pair=="Y1-Y2")] <- "Y1-Y2"
  antibody_correlation$pair2[which(antibody_correlation$pair=="O1-O2")] <-"O1-O2"
  antibody_correlation$pair2 <- factor(antibody_correlation$pair2,levels=c("Y1-Y2","O1-O2","Y-O"))
  # color <- setNames(c("#f39b7f","#f39b7f","#4dbbd5","#4dbbd5","#4dbbd5","#4dbbd5"),c("Y1-Y2","O1-O2","Y1-O1","Y1-O2","Y2-O1","Y2-O2"))
  color <- setNames(c("#f39b7f","#4dbbd5","grey"),c("Y1-Y2","O1-O2","Y-O")) 
  p_list[[antibody]] <- ggplot(antibody_correlation, aes(x = pair2, y = cor,fill=pair2)) +
    geom_boxplot(outliers =F) +
    scale_fill_manual(values=color)+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab(NULL) + ylab("Correlation")+
    ggtitle(antibody)
}
for(antibody in antibodys){
  p <- p_list[[antibody]]
  ggsave(paste0("result/Sup_figures/",antibody,"_correlation_boxplot.pdf"),p,width = 6,height = 8)
}

###WGBS
antibodys <- c("WGBS")
tissues <- c("brain","liver","testis","colon","kidney","lung","spleen","muscle","pancreas","Hip","cecum","bonemarrow","ileum","heart","thymus","stomach","skin","aorta","tongue","bladder","CB","jejunum","uterus","ovary")
bin_size <- "10kb"
for(antibody in antibodys){
  antibody_correlation <- data.frame()
  for(tissue in tissues){
    search_table <- read.csv("data/samples/all/WGBS_search_table.csv")
    search_table <- search_table[which(search_table$tissue==tissue),]
    df <- read.csv(paste0("data/samples/WGBS/",tissue,"/compress2bin/",bin_size,"_bins_all_depth.csv"))
    df <- df[which(df$total_V5 > 10),]
    df <- reshape2::dcast(df, label ~ sample, value.var = "percent")
    rownames(df) <- df$label
    
    df <-df[,-1]
    df <- df[complete.cases(df), ]
    vars <- apply(df, 1, var, na.rm = TRUE)
    df <- df[order(vars, decreasing = TRUE)[1:100000], ]

    
    search_table$age <- factor(search_table$age, levels = c("3M","24M"))
    search_table <- search_table[order(search_table$age),]
    search_table$batch <- rep(paste0("batch",c(1:(nrow(search_table)/2))),2)
    search_table <- search_table[which(search_table$batch %in% c("batch1","batch2")),]
    
    search_table$condition <- paste0(search_table$age,"_",search_table$batch)
    search_table$condition[which(search_table$condition=="3M_batch1")] <- "Y1"
    search_table$condition[which(search_table$condition=="3M_batch2")] <- "Y2"
    search_table$condition[which(search_table$condition=="24M_batch1")] <- "O1"
    search_table$condition[which(search_table$condition=="24M_batch2")] <- "O2"
    pair_mat <- combn(search_table$condition, 2)
    pair_df  <- as.data.frame(t(pair_mat))
    
    
    tissue_correaltion <- data.frame()
    for(i in 1:nrow(pair_df)){
      condition1 <- pair_df$V1[i]
      condition2 <- pair_df$V2[i]
      t_df <- df[,as.character(c(search_table$sample_name[which(search_table$condition==condition1)],search_table$sample_name[which(search_table$condition==condition2)]))]
      cor <- cor.test(t_df[,1],t_df[,2])
      t_tissue_correlation <- data.frame(tissue=tissue_label_change(tissue),cor=as.numeric(cor$estimate[1]),pair=paste0(condition1,"-",condition2))
      tissue_correaltion <- rbind(tissue_correaltion,t_tissue_correlation)
    }
    antibody_correlation <- rbind(antibody_correlation,tissue_correaltion)
  }
  antibody_correlation$pair <- factor(antibody_correlation$pair,levels=c("Y1-Y2","O1-O2","Y1-O1","Y2-O2","Y1-O2","Y2-O1"))
  antibody_correlation$pair2 <- "Y-O"
  antibody_correlation$pair2[which(antibody_correlation$pair=="Y1-Y2")] <- "Y1-Y2"
  antibody_correlation$pair2[which(antibody_correlation$pair=="O1-O2")] <-"O1-O2"
  antibody_correlation$pair2 <- factor(antibody_correlation$pair2,levels=c("Y1-Y2","O1-O2","Y-O"))
  # color <- setNames(c("#f39b7f","#f39b7f","#4dbbd5","#4dbbd5","#4dbbd5","#4dbbd5"),c("Y1-Y2","O1-O2","Y1-O1","Y1-O2","Y2-O1","Y2-O2"))
  color <- setNames(c("#f39b7f","#4dbbd5","grey"),c("Y1-Y2","O1-O2","Y-O")) 
  p_list[[antibody]] <- ggplot(antibody_correlation, aes(x = pair2, y = cor,fill=pair2)) +
    geom_boxplot(outliers =F) +
    scale_fill_manual(values=color)+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab(NULL) + ylab("Correlation")+
    ggtitle(antibody)
}
for(antibody in antibodys){
  p <- p_list[[antibody]]
  ggsave(paste0("result/Sup_figures/",antibody,"_correlation_boxplot.pdf"),p,width = 6,height = 8)
}
