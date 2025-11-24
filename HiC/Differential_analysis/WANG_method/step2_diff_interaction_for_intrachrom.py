import numpy as np
import math
from scipy.stats import binom
import os  
import gzip  
import pandas as pd  
import argparse 

parser = argparse.ArgumentParser(description='Process tissue and resolution.')  
parser.add_argument('-t', '--tissue', type=str, required=True, help='Type of tissue')  
parser.add_argument('-r', '--resolution', type=str, required=True, help='Resolution value')  
args = parser.parse_args()  
tissue = args.tissue  
resolution = args.resolution  
print(f'Tissue: {tissue}')  
print(f'Resolution: {resolution}')  

local = 200
minc = 10
search_table = pd.read_csv("/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_search_table.csv") 
search_table = search_table[search_table['tissue'] == tissue]  

young_sample = search_table[search_table['age'] == "3M"] 
young_sample = young_sample['sample_name'].tolist()
old_sample = search_table[search_table['age'] == "24M"] 
old_sample = old_sample['sample_name'].tolist()
data_path = "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/"+tissue+"/dense_matrix/"
mean_path = "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/"+tissue+"/mean_variance/"

for i in range(1,len(young_sample)+1):
    file1 = np.loadtxt(mean_path + young_sample[i-1] + "_" + resolution + ".mean_variance_number")
    lens1 = file1.shape[0]
    if i == 1:
        means = np.zeros([lens1,4])
        means[:,0] = file1[0:lens1,0]
    else:
        means[:,1] = file1[0:lens1,0]

for i in range(1,len(old_sample)+1):
    file1 = np.loadtxt(mean_path + old_sample[i-1] + "_" + resolution + ".mean_variance_number")
    lens2 = file1.shape[0]
    if lens1!=lens2:
        print('Error1 !!!!!\n')
    means[:,1+i] = file1[0:lens1,0]

re = np.zeros([int(3e7),6])
ids = 0
chrs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X']

for i in range(1,21):
    print("chr"+chrs[i-1]+" begin")
    data1 = np.loadtxt(data_path+young_sample[0]+"/raw/"+young_sample[0]+"_"+resolution+"_chr"+chrs[i-1]+"_dense.matrix")  
    data2 = np.loadtxt(data_path+young_sample[1]+"/raw/"+young_sample[1]+"_"+resolution+"_chr"+chrs[i-1]+"_dense.matrix")  
    data3 = np.loadtxt(data_path+old_sample[0]+"/raw/"+old_sample[0]+"_"+resolution+"_chr"+chrs[i-1]+"_dense.matrix")  
    data4 = np.loadtxt(data_path+old_sample[1]+"/raw/"+old_sample[1]+"_"+resolution+"_chr"+chrs[i-1]+"_dense.matrix")  
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
result_path="/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/"+tissue+"/dynamics/"
os.makedirs(result_path, exist_ok=True)  
np.savetxt(result_path+tissue+'_'+resolution+'.dynamics',re2,fmt='%.2f',delimiter='\t')