bash new_data_preprocess $tissue
#pairsqc
bash pairsqc.sh $tissue
Rscript pairsqc_log10prob_plot.R

#format transfer
bash hicpro2juicer_merge.sh $tissue
bash build_HiCpro_matrix.sh $tissue
bash sparse2dense.sh $tissue
bash distance_normalize.sh $tissue

#compartment
bash compartment/homer_compartment.sh $tissue
Rscript quality_control/homer_pca.R

#TAD
bash TAD/domaincaller.sh $tissue
bash TAD/arrowhead.sh $tissue
bash TAD/insulation_score.sh $tissue
Rscript TAD/insulation_score_tad_transfer.R

#loop
bash loop/mustache_call_loop.sh $tissue
bash loop/HiCCUP_call_loop.sh $tissue

#with Histone
Rscript with_histone_modification/histone_change_in_compartment.R 
