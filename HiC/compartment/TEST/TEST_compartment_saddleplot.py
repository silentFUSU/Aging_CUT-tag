# import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, subprocess
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from saddle_plot import saddleplot
# Import python package for working with cooler files and tools for analysis
import cooler
import cooltools.lib.plotting
import cooltools
import bioframe
import h5py 
import os  
import argparse  
parser = argparse.ArgumentParser(description='Process some integers.')  
parser.add_argument('-t', '--tissue', type=str, required=True, help='The type of tissue, e.g., kidney')  
parser.add_argument('-r', '--resolution', type=str, required=True, help='The resolution, e.g., 50000')
args = parser.parse_args()  
tissue = args.tissue  
resolution = args.resolution 
mm10_genome = bioframe.load_fasta('/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.fa')
search_table = pd.read_csv("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv") 
search_table = search_table[search_table['tissue'] == tissue]  
sample_names = search_table['sample_name'].tolist()  
os.makedirs("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/HiC/"+tissue, exist_ok=True)
os.makedirs("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/HiC/"+tissue+"/compartment/", exist_ok=True)
for sample in sample_names:
        clr = cooler.Cooler('/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/'+tissue+'/cool/'+sample+'_'+resolution+'.cool')
        with clr.open('r+') as f:
                if "SCALE_mult" in f["bins"]:
                        print(sample+" SCALE_mult exist")
                else:
                        f["bins"].create_dataset("SCALE_mult", data=1./f["bins/SCALE"][:], compression="gzip", compression_opts=6)
                
        bins = clr.bins()[:]
        gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], mm10_genome)
        view_df = pd.DataFrame({'chrom': clr.chromnames,
                                'start': 0,
                                'end': clr.chromsizes.values,
                                'name': clr.chromnames}
                        )
        # cis_eigs = cooltools.eigs_cis(
        #                         clr,
        #                         gc_cov,
        #                         view_df=view_df,
        #                         n_eigs=3,
        #                         clr_weight_name='SCALE_mult')
        # eigenvector_track = cis_eigs[1][['chrom','start','end','E1']]
        # eigenvector_track = eigenvector_track.dropna(subset=['E1'])
        
        eigenvector_track = pd.read_csv('/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/'+tissue+'/compartment/homer_compartment/PC1/'+sample+'_'+resolution+'.PC1.bedGraph', sep='\t', header=None, names=['chrom', 'start', 'end', 'E1'], skiprows=1) 
        valid_chromosomes = [f'chr{i}' for i in range(1, 20)] + ['chrX'] 
        eigenvector_track = eigenvector_track[eigenvector_track['chrom'].isin(valid_chromosomes)] 
        cvd = cooltools.expected_cis(
                clr=clr,
                view_df=view_df,
                clr_weight_name='SCALE_mult'
        )
        Q_LO = 0.025 # ignore 2.5% of genomic bins with the lowest E1 values
        Q_HI = 0.975 # ignore 2.5% of genomic bins with the highest E1 values
        N_GROUPS = 38 # divide remaining 95% of the genome into 38 equisized groups, 2.5% each
        interaction_sum, interaction_count =  cooltools.saddle(
                clr,
                cvd,
                eigenvector_track,
                'cis',
                clr_weight_name='SCALE_mult',
                n_bins=N_GROUPS,
                qrange=(Q_LO,Q_HI),
                view_df=view_df
        )
        grid = saddleplot(eigenvector_track,
                interaction_sum/interaction_count,
                N_GROUPS,
                qrange=(Q_LO,Q_HI),
                cbar_kws={'label':'average observed/expected contact frequency'}
                );
        plt.savefig("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/HiC/"+tissue+"/compartment/cooler_compartment_saddle_"+sample+"_homer.png", dpi=300)
        plt.close()
        interaction_sum_count = pd.DataFrame(interaction_sum/interaction_count)
        os.makedirs("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/"+tissue+"/compartment/cooler/", exist_ok=True)  
        interaction_sum_count.to_csv("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/"+tissue+"/compartment/cooler/"+sample+"interaction_sum_count_homer.csv",index = False)
# f, ax = plt.subplots(figsize=(15, 10))  

# norm = LogNorm(vmin=0.1,vmax=1)  

# im = ax.matshow(  
#     clr.matrix(balance=False)[:],  
#     norm=norm,  
#     cmap='fall'  
# )  
# plt.axis([0, 500, 500, 0])  

# divider = make_axes_locatable(ax)  
# cax = divider.append_axes("right", size="5%", pad=0.1)  
# plt.colorbar(im, cax=cax, label='corrected frequencies')  

# ax.set_ylabel('chr2:0-50Mb')  
# ax.xaxis.set_visible(False)  

# ax1 = divider.append_axes("top", size="20%", pad=0.25, sharex=ax)  
# weights = clr.bins()[:]['SCALE_mult'].values  
# ax1.plot([0, 500], [0, 0], 'k', lw=0.25)  
# ax1.plot(eigenvector_track['E1'].values, label='E1')  

# ax1.set_ylabel('E1')  
# ax1.set_xticks([])  

# for i in np.where(np.diff((cis_eigs[1]['E1'] > 0).astype(int)))[0]:  
#     ax.plot([0, 500], [i, i], 'k', lw=0.5)  
#     ax.plot([i, i], [0, 500], 'k', lw=0.5)  

# plt.savefig("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/HiC/lung/compartment/cooler_compartment_test.png", dpi=300)
# plt.close()

#To create saddles in cis with saddle
