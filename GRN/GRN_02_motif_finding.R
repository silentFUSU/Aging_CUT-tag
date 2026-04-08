#motif finding for each enhancer

library(motifmatchr)
library(GenomicRanges)
library(TFBSTools)
library(universalmotif)
library(edgeR)
library(Matrix)

peak=read.table('/path/to/all_enhancer.bed',header = F, sep = '\t')
#This is bed file of enhancers
rownames(peak)=paste0(peak$V1,'_',peak$V2)
peak_r=GRanges(seqnames = peak$V1,ranges = IRanges(start = peak$V2,end=peak$V3-1))

#load motif profile (.meme file)
motif_list <- read_meme("path/to/motif.meme")
pwm_list <- lapply(motif_list, function(m) {
  PWMatrix(
    ID = m@name,
    name = m@altname,
    profileMatrix = m@motif
  )})
for (i in c(1:length(pwm_list))){
  names(pwm_list)[i]=pwm_list[[i]]@name
  colnames(pwm_list[[i]]@profileMatrix)=NULL
}

#motif finding for all motifs-enhancers
motif=names(pwm_list)
match_result=list()
for(i in motif){
  motif_ix <- matchMotifs(pwm_list[[i]], peak_r, genome = "mm10",out='scores')
  match_result[[i]]=motif_ix
}
#binary result
motif_peak_bi <- do.call(cbind, lapply(names(match_result), function(motif_name) {
  motif_data <- motifMatches(match_result[[motif_name]])
  colnames(motif_data) <- motif_name
  return(motif_data)
}))
rownames(motif_peak_bi)=rownames(peak)
save(motif_peak_bi,file = '/path/to/motif/finding/result.rdata')
