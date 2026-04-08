import numpy as np
import gzip  
import os  
import pandas as pd  
import argparse 

def combine_mean_var(u1,u2,v1,v2,n1,n2):
    u3 = 1.0/(n1+n2)*(n1*u1+n2*u2)
    v3 = 1.0/(n1+n2)*(n1*(v1+u1**2-u3**2)+n2*(v2+u2**2-u3**2))
    return u3,v3

parser = argparse.ArgumentParser(description='Process tissue and resolution.')  
parser.add_argument('-t', '--tissue', type=str, required=True, help='Type of tissue')  
parser.add_argument('-r', '--resolution', type=str, required=True, help='Resolution value')  
args = parser.parse_args()  
tissue = args.tissue  
resolution = args.resolution  
print(f'Tissue: {tissue}')  
print(f'Resolution: {resolution}')  


data_path = "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/"+tissue+"/dense_matrix/"
result_path = "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/"+tissue+"/mean_variance/"
os.makedirs(result_path, exist_ok=True)  
search_table = pd.read_csv("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv") 
search_table = search_table[search_table['tissue'] == tissue]  
sample_names = search_table['sample_name'].tolist()  
for sample in sample_names: 
    print(sample+" begin")
    for i in range(1,20):
        data = np.loadtxt(data_path+sample+"/raw/"+sample+"_"+resolution+"_chr"+str(i)+"_dense.matrix")  
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

    data = np.loadtxt(data_path+sample+"/raw/"+sample+"_"+resolution+"_chrX_dense.matrix")  
    zero_row = np.sum(data,axis = 0)==0
    lens = data.shape[0]
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

    # data = np.loadtxt(data_path+sample+"/raw/"+sample+"_"+resolution+"_chrY_dense.matrix")  
    # zero_row = np.sum(data,axis = 0)==0
    # lens = data.shape[0]
    # sign = np.zeros([lens,lens])
    # sign[zero_row==1,:] = 1
    # sign[:,zero_row==1] = 1
    # for j in range(0,lens):
    #     diag1 = np.diag(data,k=j)
    #     diag2 = np.diag(sign,k=j)
    #     diagu = diag1[diag2==0]
    #     if diagu.shape[0]!=0:
    #         u2 = np.mean(diagu)
    #         v2 = np.mean(diagu)
    #         n2 = diagu.shape[0]
    #         if m_v[j,2] == 0:
    #             m_v[j,0:3] = np.array([u2,v2,n2])
    #         else:
    #             u3,v3 = combine_mean_var(m_v[j,0],u2,m_v[j,1],v2,m_v[j,2],n2)
    #             m_v[j,0:3] = np.array([u3,v3,m_v[j,2]+n2])  

    np.savetxt(result_path+sample+'_'+resolution+'.mean_variance_number',m_v,fmt='%.3f',delimiter='\t')
    print(sample+" done")

print("all done")
