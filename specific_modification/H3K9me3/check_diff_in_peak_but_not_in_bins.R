library(data.table)
tissue <- "ileum"
df <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_young_old_merge-W1000-G3000-E100_diff_after_remove_batch_effect.csv"))
bin <- read.csv(paste0("data/samples/",tissue,"/H3K9me3/H3K9me3_10kb_bins_diff_after_remove_batch_effect.csv"))
bin <- bin[which(bin$Significant=="Down"),]
bin <- bin[,c("Chr","Start","End")]
colnames(bin) <- c("V1","V2","V3")
df <- df[which(df$Significant=="Down"),]
peaks <- df[,c("Chr","Start","End")]
peaks <- as.data.table(peaks)
setDT(peaks)
setkey(peaks,"Chr","Start","End")

bin <- as.data.table(bin)
setDT(bin)
setkey(bin,"V1","V2","V3")
overlaps <- foverlaps(bin,peaks, type = "any", nomatch = 0L)
overlaps$Geneid <- paste0(overlaps$V1,":",overlaps$Start,"-",overlaps$End)

df <- df[-which(df$Geneid %in% overlaps$Geneid),]
