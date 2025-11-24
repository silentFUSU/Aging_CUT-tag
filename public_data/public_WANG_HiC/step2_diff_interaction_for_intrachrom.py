## DESCRIPTIONS of THIS CODE (only used for intrachromosome)
# author: X.Z.
# reference Dixon et al. 2012
# (STEP 2)This code was used to perform differential interaction analysis and calculate P values, and P values of background. (See methods)
#input files:
#		1.the list of raw matrix (not be normalized) files (N*N) from chr1:chrX; two replicates of two states (4 lists)
#		2.mean_variance_number which calculated in last step. we need four files for two replicates of two states. (4 files)
#output files:
#		1.the list of P values for one foreground and two background


#note: we have mark the positions where we need to set again using "##**##".


sample1 = 'G'##**##state 1
sample2 = 'DS'##**##state 2
path1 = '/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/WANG_cellular_aging_HiC/mean_variance/' ##**##basic path
resolution="20000"
local = 300#(6M/20k) ##**##the max distance(/bin) to calculate . This is depend on matrix resolution and sequencing depth (40M/200k)
minc = 10#(mini counts) ##**## sum of count for one interaction less than this num will be ignored
import numpy as np
import math
from scipy.stats import binom
#means
for i in range(1,3):
    file1 = np.loadtxt(path1+sample1+str(i)+'_'+resolution+'.mean_variance_number')##**##two replicates mean_variance_number files
    lens1 = file1.shape[0]
    if i == 1:
        means = np.zeros([lens1,4])
        means[:,0] = file1[0:lens1,0]
    else:
        means[:,1] = file1[0:lens1,0]
for i in range(1,3):
    file1 = np.loadtxt(path1+sample2+str(i)+'_'+resolution+'.mean_variance_number')##**##two replicates mean_variance_number files
    lens2 = file1.shape[0]
    if lens1!=lens2:
        print('Error1 !!!!!\n')
    means[:,1+i] = file1[0:lens1,0]
#dynamic interactions
re = np.zeros([int(3e7), 6])  
ids = 0
path1="/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/WANG_cellular_aging_HiC/dense_matrix/"
for i in range(1,23):
    data1 = np.loadtxt(path1+sample1+'1/raw/'+sample1+'1_'+resolution+'_chr'+str(i)+'_dense.matrix')##**## replicate1 matrix of state 1
    data2 = np.loadtxt(path1+sample1+'2/raw/'+sample1+'2_'+resolution+'_chr'+str(i)+'_dense.matrix')##**## replicate2 matrix of state 1
    data3 = np.loadtxt(path1+sample2+'1/raw/'+sample2+'1_'+resolution+'_chr'+str(i)+'_dense.matrix')##**## replicate1 matrix of state 2
    data4 = np.loadtxt(path1+sample2+'2/raw/'+sample2+'2_'+resolution+'_chr'+str(i)+'_dense.matrix')##**## replicate2 matrix of state 2
    lens1 = data1.shape[0]
    zrow = np.zeros([lens1,4])
    zrow[:,0] = (np.sum(data1,axis = 1)==0)
    zrow[:,1] = (np.sum(data2,axis = 1)==0)
    zrow[:,2] = (np.sum(data3,axis = 1)==0)
    zrow[:,3] = (np.sum(data4,axis = 1)==0)
    for j in range(0,lens1):
        if sum(zrow[j,0:4])==0:
            s = max(0,j-local)
            o = j
            for k in range(s,o+1):
                val = [data1[j,k],data2[j,k],data3[j,k],data4[j,k]]
                if sum(zrow[k,0:4])==0 and sum(val)>=minc:
                    p1 = sum(means[j-k,0:2])*1.0/sum(means[j-k,0:4])
                    p2 = sum(means[j-k,[0,2]])*1.0/sum(means[j-k,0:4])
                    p3 = sum(means[j-k,[0,3]])*1.0/sum(means[j-k,0:4])
                    trials = math.ceil(sum(val))
                    if sum(val[0:2])-1>=p1*trials:
                        pval1 = math.log10(1-binom.cdf(sum(val[0:2])-1,trials,p1)+1e-5)
                    else:
                        pval1 = -1*math.log10(binom.cdf(sum(val[0:2]),trials,p1)+1e-5)
                    if val[0]+val[2]-1>=p2*trials:
                        pval2 = math.log10(1-binom.cdf(val[0]+val[2]-1,trials,p2)+1e-5)
                    else:
                        pval2 = -1*math.log10(binom.cdf(val[0]+val[2],trials,p2)+1e-5)
                    if val[0]+val[3]-1>=p3*trials:
                        pval3 = math.log10(1-binom.cdf(val[0]+val[3]-1,trials,p3)+1e-5)
                    else:
                        pval3 = -1*math.log10(binom.cdf(val[0]+val[3],trials,p3)+1e-5)
                    re[ids,:] = [i,j,k,round(pval1,2),round(pval2,2),round(pval3,2)]
                    ids = ids+1

#chrX
for i in range(23,24):
    data1 = np.loadtxt(path1+sample1+'1/raw/'+sample1+'1_'+resolution+'_chrX'+'_dense.matrix')
    data2 = np.loadtxt(path1+sample1+'2/raw/'+sample1+'2_'+resolution+'_chrX'+'_dense.matrix')
    data3 = np.loadtxt(path1+sample2+'1/raw/'+sample2+'1_'+resolution+'_chrX'+'_dense.matrix')
    data4 = np.loadtxt(path1+sample2+'2/raw/'+sample2+'2_'+resolution+'_chrX'+'_dense.matrix')
    lens1 = data1.shape[0]
    zrow = np.zeros([lens1,4])
    zrow[:,0] = (np.sum(data1,axis = 1)==0)
    zrow[:,1] = (np.sum(data2,axis = 1)==0)
    zrow[:,2] = (np.sum(data3,axis = 1)==0)
    zrow[:,3] = (np.sum(data4,axis = 1)==0)
    for j in range(0,lens1):
        if sum(zrow[j,0:4])==0:
            s = max(0,j-local)
            o = j
            for k in range(s,o+1):
                val = [data1[j,k],data2[j,k],data3[j,k],data4[j,k]]
                if sum(zrow[k,0:4])==0 and sum(val)>=minc:
                    p1 = sum(means[j-k,0:2])*1.0/sum(means[j-k,0:4])
                    p2 = sum(means[j-k,[0,2]])*1.0/sum(means[j-k,0:4])
                    p3 = sum(means[j-k,[0,3]])*1.0/sum(means[j-k,0:4])
                    trials = math.ceil(sum(val))
                    if sum(val[0:2])-1>=p1*trials:
                        pval1 = math.log10(1-binom.cdf(sum(val[0:2])-1,trials,p1)+1e-5)
                    else:
                        pval1 = -1*math.log10(binom.cdf(sum(val[0:2]),trials,p1)+1e-5)
                    if val[0]+val[2]-1>=p2*trials:
                        pval2 = math.log10(1-binom.cdf(val[0]+val[2]-1,trials,p2)+1e-5)
                    else:
                        pval2 = -1*math.log10(binom.cdf(val[0]+val[2],trials,p2)+1e-5)
                    if val[0]+val[3]-1>=p3*trials:
                        pval3 = math.log10(1-binom.cdf(val[0]+val[3]-1,trials,p3)+1e-5)
                    else:
                        pval3 = -1*math.log10(binom.cdf(val[0]+val[3],trials,p3)+1e-5)
                    re[ids,:] = [i,j,k,round(pval1,2),round(pval2,2),round(pval3,2)]
                    ids = ids+1
re2 = re[0:ids,:]
np.savetxt('/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/WANG_cellular_aging_HiC/dynamics/'+sample1+'_'+sample2+'_'+resolution+'.dynamics',re2,fmt='%.2f',delimiter='\t')
                        
                    
                    
                    
                