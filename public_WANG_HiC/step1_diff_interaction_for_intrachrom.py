## DESCRIPTIONS of THIS CODE(only used for intrachromosome)
# author: X.Z.
# reference Dixon et al. 2012
# (STEP 1)This code was used to calculate mean and variance of contacts across different genomic distances, which is used for calculate the differential interactions
#input files:
#		1.the list of raw matrix (not be normalized) files (N*N) from chr1:chrX; matrices of each replicate need to be calculated separately;
#output files:
#		1.the list of mean and variance of contact across different genomic distances at whole-genome level


#note: if you have 2 replicates for two states, you need to have 4 output files to do STEP 2.
#note: we have mark the positions where we need to set again using "##**##".


sample = 'G2'

def combine_mean_var(u1,u2,v1,v2,n1,n2):
    u3 = 1.0/(n1+n2)*(n1*u1+n2*u2)
    v3 = 1.0/(n1+n2)*(n1*(v1+u1**2-u3**2)+n2*(v2+u2**2-u3**2))
    return u3,v3
import numpy as np
#####deal with chrom 1:22
data_path="/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/WANG_cellular_aging_HiC/dense_matrix/"+sample+"/"
resolution="200000"
for i in range(1,23):##**##chromosome number you used
    data = np.loadtxt(data_path+'/raw/'+sample+'_'+resolution+'_chr'+str(i)+'_dense.matrix') ##**##matrix path and name
    zero_row = np.sum(data,axis = 0)==0
    lens = data.shape[0]
    #build sign matrix
    sign = np.zeros([lens,lens])
    sign[zero_row==1,:] = 1
    sign[:,zero_row==1] = 1
    if i == 1:
        m_v = np.zeros([lens,3])
    for j in range(0,lens):
        diag1 = np.diag(data,k=j)
        diag2 = np.diag(sign,k=j)
        diagu = diag1[diag2==0]
        if diagu.shape[0]!=0:
            u2 = np.mean(diagu)
            v2 = np.var(diagu)
            n2 = diagu.shape[0]
            if m_v[j,2] == 0:
                m_v[j,0:3] = np.array([u2,v2,n2])
            else:
                u3,v3 = combine_mean_var(m_v[j,0],u2,m_v[j,1],v2,m_v[j,2],n2)
                m_v[j,0:3] = np.array([u3,v3,m_v[j,2]+n2])
####deal with chrom X
data = np.loadtxt(data_path+'/raw/'+sample+'_'+resolution+'_chrX_dense.matrix') ##**##if you wanted to used chrX 
zero_row = np.sum(data,axis = 0)==0
lens = data.shape[0]
#build sign matrix
sign = np.zeros([lens,lens])
sign[zero_row==1,:] = 1
sign[:,zero_row==1] = 1
for j in range(0,lens):
    diag1 = np.diag(data,k=j)
    diag2 = np.diag(sign,k=j)
    diagu = diag1[diag2==0]
    if diagu.shape[0]!=0:
        u2 = np.mean(diagu)
        v2 = np.mean(diagu)
        n2 = diagu.shape[0]
        if m_v[j,2] == 0:
            m_v[j,0:3] = np.array([u2,v2,n2])
        else:
            u3,v3 = combine_mean_var(m_v[j,0],u2,m_v[j,1],v2,m_v[j,2],n2)
            m_v[j,0:3] = np.array([u3,v3,m_v[j,2]+n2])   
np.savetxt('/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/WANG_cellular_aging_HiC/mean_variance/'+sample+'_'+resolution+'.mean_variance_number',m_v,fmt='%.3f',delimiter='\t')
