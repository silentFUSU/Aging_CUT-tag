# -*- coding: utf-8 -*-
"""
@author: X.Z.
transform BedGraph to bin levels
"""


import pandas as pd
import numpy as np
import sys,getopt,math
res = 1e4

opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["inputfile=","outputfile="])
for opt, arg in opts:
    if opt == '-h':
        print ('python bedtobins.py -i [inputfile] -o [outputfile] ')
        sys.exit()
    elif opt in ("-i", "--inputfile"):
        path1 = arg
    elif opt in ("-o","--outputfile"):
        path2 = arg




chrom_size = pd.read_csv('/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.normal.chrom.sizes',sep = '\t',header=None)
signal = pd.read_csv(path1,sep = '\t',header = None)
signal.columns = ['chrom','start','end','fc']


merge = np.zeros((0,4))
for i in range (0,chrom_size.shape[0]):
    chrom = chrom_size.iloc[i,0]
    print(chrom)
    chrom_n = chrom_size.iloc[i,1]
    nums = int(math.ceil(chrom_n/res))
    chrom_bins = np.zeros((nums,3))
    chrom_bins[0:nums,0] = range(1,int(res*nums+1),int(res))
    chrom_bins[0:nums-1,1] = range(int(res),int(res*nums),int(res))
    chrom_bins[nums-1,1] = chrom_n
    signal_chrom0 = signal.loc[signal.iloc[:,0]==chrom,:]
    signal_chrom = signal_chrom0.astype({'start': 'int64','end': 'int64','fc': 'float64'})
    signal_chrom.sort_values(by='start')
    len1 = signal_chrom.shape[0]
    len2 = chrom_bins.shape[0]
    x = 0
    y = 0
    while x<len1 and y<len2:
        if signal_chrom.iloc[x,1]+1>=chrom_bins[y,0] and signal_chrom.iloc[x,2]<=chrom_bins[y,1]:
            chrom_bins[y,2] += signal_chrom.iloc[x,3]*(signal_chrom.iloc[x,2]-signal_chrom.iloc[x,1])
            x += 1
        elif signal_chrom.iloc[x,1]+1>chrom_bins[y,1]:
            y += 1
        elif signal_chrom.iloc[x,2]<chrom_bins[y,0]:
            x += 1
        elif signal_chrom.iloc[x,1]+1>=chrom_bins[y,0] and signal_chrom.iloc[x,2]>chrom_bins[y,1]:
            chrom_bins[y,2] += signal_chrom.iloc[x,3]*(chrom_bins[y,1]-signal_chrom.iloc[x,1]+1)
            y += 1
        elif signal_chrom.iloc[x,1]+1<chrom_bins[y,0] and signal_chrom.iloc[x,2]<=chrom_bins[y,1]:
            chrom_bins[y,2] += signal_chrom.iloc[x,3]*(signal_chrom.iloc[x,2]-chrom_bins[y,0])
            x += 1
        elif signal_chrom.iloc[x,1]+1<chrom_bins[y,0] and signal_chrom.iloc[x,2]>chrom_bins[y,1]:
            chrom_bins[y,2] += signal_chrom.iloc[x,3]*(chrom_bins[y,1]-chrom_bins[y,0]+1)
            y += 1
    chrom2 = chrom.split('r')[1]
    if chrom2=='X':
        chrom2 = 23
    else:
        chrom2 = int(chrom2)
    chroms_part = np.ones((chrom_bins.shape[0],1))*chrom2
    chroms_out = np.concatenate((chroms_part,chrom_bins),axis =1)
    merge = np.concatenate((merge,chroms_out),axis = 0)
merge[:,3] = merge[:,3]/(merge[:,2]-merge[:,1]+1)
np.savetxt(path2,merge,fmt='%d %d %d %.3f',delimiter = '\t')
