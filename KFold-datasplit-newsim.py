# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 11:02:27 2019

@author: Hp
"""
import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns
from tqdm import tqdm
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RepeatedKFold
import os 
#读取#读取miRNA-Drug关系矩阵，进行5-fold交叉验证前的数据集划分
assoc_path = u'../data/miRNA-Drug association/'
#pro_assoc_array = pd.read_csv(assoc_path+'Pro-mTD/mir-drug-association.csv')
ncdr_assoc_array = pd.read_csv(assoc_path+'ncDR/mir-drug-array.csv')
kf_assoc_path = u'../data/KFold/5-fold/'

K=5
rkf = RepeatedKFold(n_splits=K, n_repeats=500)
i=1
for train, test in rkf.split(ncdr_assoc_array):
    train_kf = ncdr_assoc_array.loc[train,]
    train_kf.to_csv(kf_assoc_path+'ncDR/train/'+str(i)+'.csv',index=None)
    i=i+1

#miRNA-name和Drug-name
drug_name = pd.read_csv('../data/Drug/DrugName.txt',sep='\t',header=None)
new_col = ['Drug']
drug_name.columns= new_col
mir_name = pd.read_csv('../data/miRNA/miRName.txt',sep='\t',header=None)
new_col = ['miRNA']
mir_name.columns= new_col
banch = 50000
#转为关系矩阵
def array_matrix(df):
    df = pd.DataFrame(df)
    temp = list(zip(df.iloc[:,1],df.iloc[:,0],df.iloc[:,2]))
    G = nx.Graph()
    G.add_weighted_edges_from(temp)
    df_matrix = nx.to_pandas_adjacency(G)
    result = df_matrix.loc[mir_name.miRNA.tolist(),drug_name.Drug.tolist()]
    return result

#生成新的相似矩阵
def new_simliar(df):
    df = pd.DataFrame(df)
    m,n = df.shape
    inds = np.array(df.index)
    df1_index = np.tile(inds.reshape(-1,1),m).ravel()
    df2_index = np.tile(inds.reshape(1,-1),m).ravel()
    df1_value = df.values.repeat(m,axis=0)
    df2_value = np.tile(df.values,(m, 1))
    ecope = int((m*m)/banch)
    values = np.zeros(m*m)
    for j in range(ecope):
        start = j * banch
        end = (j+1) * banch
        values[start:end] = (df1_value[start:end] * df2_value[start:end]).sum(axis=1)
    values[ecope * banch: m * m] = (df1_value[ecope * banch: m * m] * df2_value[ecope * banch: m * m]).sum(axis=1) 
    G = nx.Graph()
    G.add_weighted_edges_from(list(zip(df1_index,df2_index,values)))
    new_sim = nx.to_pandas_adjacency(G)
    return new_sim

    
    

for i in tqdm(range(1,500*K+1)):
    train_array = pd.read_csv(kf_assoc_path+'ncDR/train/'+str(i)+'.csv')
    #生成训练集矩阵
    assoc_data_1 = pd.merge(train_array,mir_name.miRNA,how='inner',on=['miRNA'])
    assoc_data_2 = pd.merge(assoc_data_1,drug_name.Drug,how='inner',on=['Drug'])
    assoc_data_3 = pd.merge(assoc_data_2,mir_name.miRNA,how='outer',on=['miRNA'])
    assoc_data_3['Drug']=assoc_data_3['Drug'].fillna(value='Cisplatin')
    assoc_data_3['Control']=assoc_data_3['Control'].fillna(value=0)
    assoc_data_4 = pd.merge(assoc_data_3,drug_name.Drug,how='outer',on=['Drug'])
    assoc_data_4['miRNA']=assoc_data_4['miRNA'].fillna(value='hsa-let-7')
    assoc_data_4['Control']=assoc_data_4['Control'].fillna(value=0)
    del assoc_data_1,assoc_data_2,assoc_data_3,train_array
    mir_drug_matrix = array_matrix(assoc_data_4)
    mir_drug_matrix.to_csv(kf_assoc_path+'ncDR/mir-drug-matrix/'+str(i)+'.txt',sep='\t',index=None,header=False)
    #新相似矩阵
    #new-mir-sim
    new_mir_sim = new_simliar(mir_drug_matrix)
    new_mir_sim_matrix = new_mir_sim.loc[mir_name.miRNA.tolist(),mir_name.miRNA.tolist()]
    new_mir_sim_matrix.to_csv(kf_assoc_path+'ncDR/new-mir-sim/'+str(i)+'.txt',sep='\t',index=None,header=False)
    #new-drug-sim
    new_drug_sim = new_simliar(mir_drug_matrix.T)
    new_drug_sim_matrix = new_drug_sim.loc[drug_name.Drug.tolist(),drug_name.Drug.tolist()]
    new_drug_sim_matrix.to_csv(kf_assoc_path+'ncDR/new-drug-sim/'+str(i)+'.txt',sep='\t',index=None,header=False)
    

