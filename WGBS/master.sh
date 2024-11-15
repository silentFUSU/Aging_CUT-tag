bash new_data_preprocess.sh

#quality control
/usr/local/lib64/R/bin/Rscript quality_control/CpG_overview.R

#compress to bin
bash master_CpG_compress_to_bin.sh 

#differential analysis
bash DMR/master_DSS_DMR.sh
/usr/local/lib64/R/bin/Rscript DMR_split.R
bash DMR/master_delta_methylation_1kb_bin.sh