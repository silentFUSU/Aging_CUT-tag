rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-pc-linux-gnu-library/4.2/","/usr/local/lib64/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/")
set.seed(1)
library(exomePeak2)
library(tidyr)
library(stringr)
library(dplyr)
library(GenomicFeatures)
library(TxDb.Mfascicularis.UCSC.macFas5.ncbiRefSeq)
table <- read.csv("data/public_data/m6a_monkey_CRA005942/CRA005942.csv",row.names = 1)
GENE_ANNO_GTF = "~/ref_data/Macaca_fascicularis.Macaca_fascicularis_5.0.102.gff3"
data_path = "~/projects/Aging_CUT_Tag/data/public_data/m6a_monkey_CRA005942/CRA005942/"

input_data <- table[which(str_detect(table$Run.title,"INPUT")),]
heart_input_data <- input_data[which(str_detect(input_data$Run.title,"H-") & str_detect(input_data$Run.title,"M")),]
muscle_input_data <- input_data[which(str_detect(input_data$Run.title,"M-") & !str_detect(input_data$Run.title,"F")),]

IP_data <- table[which(str_detect(table$Run.title,"IP")),]
heart_IP_data <- IP_data[which(str_detect(IP_data$Run.title,"H-") & !str_detect(IP_data$Run.title,"F")),]

heart_young_IP <- heart_IP_data$Accession[which(str_detect(heart_IP_data$Run.title,"Y"))]
heart_old_IP <- heart_IP_data$Accession[which(str_detect(heart_IP_data$Run.title,"O"))]
bam_IP <- c(paste0(data_path,heart_young_IP,"/",heart_young_IP,".bam"))
bam_IP_treated <- c(paste0(data_path,heart_old_IP,"/",heart_old_IP,".bam"))

heart_young_input <- heart_input_data$Accession[which(str_detect(heart_input_data$Run.title,"Y"))]
heart_old_input <- heart_input_data$Accession[which(str_detect(heart_input_data$Run.title,"O"))]
bam_input <- c(paste0(data_path,heart_young_input,"/",heart_young_input,".bam"))
bam_input_treated <- c(paste0(data_path,heart_old_input,"/",heart_old_input,".bam"))
txdb <- TxDb.Mfascicularis.UCSC.macFas5.ncbiRefSeq
peak_diff <- exomePeak2(bam_ip = bam_IP,
                        bam_input = bam_input ,
                        bam_ip_treated = bam_IP_treated,
                        bam_input_treated = bam_input_treated,
                        txdb=txdb,
                        genome = "BSgenome.Mfascicularis.NCBI.5.0")

