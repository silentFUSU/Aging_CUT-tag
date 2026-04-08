#aging-related TF selection

library(limma)

reg_score=read.table('path/to/pagerank/result.txt')
#This file is matrix of pagerank raw result (TF x sample)
reg_score=log(reg_score+0.000001)

#limma
all_results <- list()
all_sig <- list()
for (tissue in tissues) {
  #samples_tissue is sampleID belongs to processing tissue
  mat_tissue <- reg_score[, samples_tissue] 
  #Age is sample age information (a vector of 'young' and 'old')
  group <- factor(Age, levels = c("young", "old"))
  design <- model.matrix(~0 + group)
  colnames(design) <- c("young", "old")
  fit <- lmFit(mat_tissue, design)
  contrast <- makeContrasts(old - young, levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  results <- topTable(fit2, coef = 1, n = Inf, sort.by = "P")
  results$TF <- rownames(results)
  results$tissue <- tissue
  all_results[[tissue]] <- results
  results_sig=results[results$FDR<0.1,]
  all_sig[[tissue]]=results_sig
}
save(all_sig,file='path/to/limma/result.rdata')
