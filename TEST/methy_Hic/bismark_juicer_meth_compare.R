rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(stringr)
bismark_readname <- read.table("data/raw_data/20240804_WGBS/compare/XX315/WGBS_readnames_sorted.txt")
juicer_readname <- read.table("data/raw_data/20240804_WGBS/compare/XX315/juicer_meth_readnames_sort.txt")
juicer_unique_readname <- as.data.frame(juicer_readname[which(!juicer_readname$V1 %in% bismark_readname$V1),])
write.table(juicer_unique_readname,"data/raw_data/20240804_WGBS/compare/XX315/juicer_meth_unique_readnames_sort.txt",row.names = F,col.names = F,quote = F)

juicer_unique_methylation <- read.delim("data/raw_data/20240804_WGBS/compare/XX315/juicer_meth_unique_reads_CpG.bedGraph",skip=1,header = F)
(sum(juicer_unique_methylation$V5)/(sum(juicer_unique_methylation$V5)+sum(juicer_unique_methylation$V6)))

calculate_m_proportion <- function(cigar_string) {  
  # 提取所有的数字  
  all_numbers <- as.numeric(str_extract_all(cigar_string, "\\d+")[[1]])  
  
  # 提取M之前的数字  
  m_numbers <- as.numeric(str_extract_all(cigar_string, "\\d+(?=M)")[[1]])  
  
  # 计算总和  
  total_sum <- sum(all_numbers)  
  m_sum <- sum(m_numbers)  
  
  # 计算比例  
  proportion <- m_sum / total_sum  
  return(proportion)  
}  
bismark_read_info <- read.delim("data/raw_data/20240804_WGBS/compare/XX315/WGBS_reads_info.txt")
bismark_read_info <- bismark_read_info[order(bismark_read_info$ReadName),]
max_distance <- 0
bismark_read_info$distance <- 0
for(i in seq(1, nrow(bismark_read_info), by=2)){
  bismark_read_info[i,"distance"] <- abs(bismark_read_info[i,"Position"] - bismark_read_info[i+1,"Position"])
  bismark_read_info[i+1,"distance"] <- abs(bismark_read_info[i,"Position"] - bismark_read_info[i+1,"Position"])
  if(abs(bismark_read_info[i,"Position"] - bismark_read_info[i+1,"Position"]) > max_distance){
    max_distance <- abs(bismark_read_info[i,"Position"] - bismark_read_info[i+1,"Position"])
  }
}



df <- read.delim("data/raw_data/20240804_WGBS/compare/XX315/juicer_unique_reads_info.txt")
df$Mpercent <- sapply(df$CIGAR, calculate_m_proportion)
map_larger_90 <- df[which(df$Mpercent>0.9),]
read_name <- as.data.frame(table(map_larger_90$ReadName))
map_larger_90 <- map_larger_90[which(map_larger_90$ReadName %in% read_name$Var1[which(read_name$Freq==2)]),]

write.table(unique(map_larger_90$ReadName),"data/raw_data/20240804_WGBS/compare/XX315/map_larger_90.txt",quote = F,row.names = F,col.names = F)

df <- df[which(df$Mpercent==1),] 
read_name <- as.data.frame(table(df$ReadName))
df <- df[which(df$ReadName %in% read_name$Var1[which(read_name$Freq==2)]),]
df <- df[order(df$ReadName),]
write.table(unique(df$ReadName),"data/raw_data/20240804_WGBS/compare/XX315/all_150M.txt",quote = F,row.names = F,col.names = F)

reads_pair <- nrow(df)/2
strand_same_count=0
long_distance_count=0
dif_chr_count=0
normal_count=0

dif_chr_reads <- vector()
strand_same_reads <- vector()
long_distance_reads <- vector()
normal_reads <- vector()
for(i in seq(1, nrow(df), by=2)){
  if(df[i,"Chromosome"]!= df[i+1,"Chromosome"]){
    dif_chr_count= dif_chr_count+1
    dif_chr_reads <- c(dif_chr_reads,df[i,"ReadName"])
  }else{
    if(df[i,"Strand"] == df[i+1,"Strand"]){
      strand_same_count <- strand_same_count+1
      strand_same_reads <- c(strand_same_reads,df[i,"ReadName"])
    }else{
      if(abs(df[i,"Position"]-df[i+1,"Position"])>500){
        long_distance_count=long_distance_count+1
        long_distance_reads <- c(long_distance_reads,df[i,"ReadName"])
      }else{
        normal_count <- normal_count + 1
        normal_reads <- c(normal_reads,df[i,"ReadName"])
      }
    }
  }
}
write.table(dif_chr_reads,"data/raw_data/20240804_WGBS/compare/XX315/juicer_unique_150M_dif_chr.txt",col.names = F,quote = F,row.names = F)  
write.table(strand_same_reads,"data/raw_data/20240804_WGBS/compare/XX315/juicer_unique_150M_same_chr_same_strand.txt",col.names = F,quote = F,row.names = F)  
write.table(long_distance_reads,"data/raw_data/20240804_WGBS/compare/XX315/juicer_unique_150M_same_chr_dif_strand_larger_500bp.txt",col.names = F,quote = F,row.names = F)  
write.table(normal_reads,"data/raw_data/20240804_WGBS/compare/XX315/juicer_unique_150M_same_chr_dif_strand_smaller_500bp.txt",col.names = F,quote = F,row.names = F)  


write.table(read_name_keep,"data/raw_data/20240804_WGBS/compare/XX315/juicer_unique_150M_short_distance_same_chr_dif_strand.txt",col.names = F,quote = F,row.names = F)  

df_filter <- df[which(df$ReadName%in% read_name_keep$read_name_keep),]
  
df <- read.delim("data/raw_data/20240804_WGBS/compare/XX315/juicer_unique_reads_info.txt")
df$Mpercent <- sapply(df$CIGAR, calculate_m_proportion)
map_smaller_90 <- df[which(df$Mpercent < 0.9),]
read_name <- as.data.frame(table(map_smaller_90$ReadName))
map_smaller_90 <- map_smaller_90[which(map_smaller_90$ReadName %in% read_name$Var1[which(read_name$Freq==2)]),]
write.table(unique(map_smaller_90$ReadName),"data/raw_data/20240804_WGBS/compare/XX315/juicer_unique_Mpercent_smaller_90.txt",col.names = F,quote = F,row.names = F)  

map_larger90_smaller100 <-  df[which(df$Mpercent > 0.9 & df$Mpercent < 1),]
read_name <- as.data.frame(table(map_larger90_smaller100$ReadName))
map_larger90_smaller100 <- map_larger90_smaller100[which(map_larger90_smaller100$ReadName %in% read_name$Var1[which(read_name$Freq==2)]),]

write.table(unique(map_larger90_smaller100$ReadName),"data/raw_data/20240804_WGBS/compare/XX315/juicer_unique_Mpercent_larger90_smaller100.txt",col.names = F,quote = F,row.names = F)  

read_name <- as.data.frame(table(df$ReadName))
multi_fragment <- df[which(df$ReadName %in% read_name$Var1[which(read_name$Freq!=2)]),]
write.table(unique(multi_fragment$ReadName),"data/raw_data/20240804_WGBS/compare/XX315/juicer_unique_NOT_only_two_end.txt",col.names = F,quote = F,row.names = F)  





