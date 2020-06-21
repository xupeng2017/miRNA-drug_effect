
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 10:55:15 2019

@author: Hp
"""

import numpy as np
import pandas as pd
import networkx as nx
import os 
from tqdm import tqdm

#读取drug-drug相似数组，转为相似矩阵
drug_sim_array = pd.read_excel('../data/Drug/drug-similarity-array.xlsx',sheet_name = 'Sheet1')
temp = list(zip(drug_sim_array.iloc[:,0],drug_sim_array.iloc[:,1],drug_sim_array.iloc[:,2]))
G = nx.Graph()
G.add_weighted_edges_from(temp)
drug_sim_matrix = nx.to_pandas_adjacency(G)
drug_sim_matrix.to_csv('../data/Drug/drug-sim-matrix.csv')
drug_name = pd.DataFrame({'Drug':drug_sim_matrix.index})
drug_name.to_csv('../data/Drug/DrugName.txt',index=None,header=False,sep='\t')
for i in range(drug_sim_matrix.shape[0]):
    drug_sim_matrix.iloc[i,i] = 1
drug_sim_matrix.to_csv('../data/Drug/drug-sim-matrix.txt',index=None,header=False,sep='\t')
#读取miRNA-miRNA相似矩阵
mir_sim_matrix = pd.read_csv('../data/miRNA/miR-sim-all.csv',index_col='miRNA')
for i in range(mir_sim_matrix.shape[0]):
    mir_sim_matrix.iloc[i,i] = 1
mir_name = pd.DataFrame({'miRNA':mir_sim_matrix.index})
mir_name.to_csv('../data/miRNA/miRName.txt',index=None,header=False,sep='\t')
mir_sim_matrix.to_csv('../data/miRNA/mir-sim-matrix.txt',index=None,header=False,sep='\t')
del drug_sim_array,temp,G

assoc_path = u'../data/miRNA-Drug association/'
#读取pro-mTD_miR-drug关系网络生成关系矩阵
pro_assoc_array = pd.read_csv(assoc_path+'Pro-mTD/mir-drug-association.csv')
assoc_data_1 = pd.merge(pro_assoc_array,mir_name.miRNA,how='inner',on=['miRNA'])
assoc_data_2 = pd.merge(assoc_data_1,drug_name.Drug,how='inner',on=['Drug'])
assoc_data_3 = pd.merge(assoc_data_2,mir_name.miRNA,how='outer',on=['miRNA'])
assoc_data_3['Drug']=assoc_data_3['Drug'].fillna(value='Cisplatin')
assoc_data_3['Control']=assoc_data_3['Control'].fillna(value=0)
assoc_data_4 = pd.merge(assoc_data_3,drug_name.Drug,how='outer',on=['Drug'])
assoc_data_4['miRNA']=assoc_data_4['miRNA'].fillna(value='hsa-let-7')
assoc_data_4['Control']=assoc_data_4['Control'].fillna(value=0)
del assoc_data_1,assoc_data_2,assoc_data_3,pro_assoc_array
temp = list(zip(assoc_data_4.iloc[:,1],assoc_data_4.iloc[:,0],assoc_data_4.iloc[:,2]))
G = nx.Graph()
G.add_weighted_edges_from(temp)
pro_mir_drug_matrix = nx.to_pandas_adjacency(G)
pro_mir_drug_matrix = pro_mir_drug_matrix.loc[mir_name.miRNA.tolist(),drug_name.Drug.tolist()]
pro_mir_drug_matrix.to_csv(assoc_path+'Pro-mTD/mir-drug-matrix.csv')
pro_mir_drug_matrix.to_csv(assoc_path+'Pro-mTD/mir-drug-matrix.txt',sep='\t',index=None,header=False)
del assoc_data_4,temp,G
#生成新的相似矩阵
m,n = pro_mir_drug_matrix.shape
banch = 50000
#new-mirna
inds = np.array(pro_mir_drug_matrix.index)
df1_index = np.tile(inds.reshape(-1,1),m).ravel()
df2_index = np.tile(inds.reshape(1,-1),m).ravel()
df1_value = pro_mir_drug_matrix.values.repeat(m,axis=0)
df2_value = np.tile(pro_mir_drug_matrix.values,(m, 1))
ecope = int((m*m)/banch)
values = np.zeros(m*m)
for i in range(ecope):
    start = i * banch
    end = (i+1) * banch
    values[start:end] = (df1_value[start:end] * df2_value[start:end]).sum(axis=1)         
values[ecope * banch: m * m] = (df1_value[ecope * banch: m * m] * df2_value[ecope * banch: m * m]).sum(axis=1) 
G = nx.Graph()
G.add_weighted_edges_from(list(zip(df1_index,df2_index,values)))
mir_sim_new = nx.to_pandas_adjacency(G)
mir_sim_new_matrix = mir_sim_new.loc[mir_name.miRNA.tolist(),mir_name.miRNA.tolist()]
mir_sim_new_matrix.to_csv('../data/new similarity/Pro-mTD/new-mir-sim.txt',header=None,index=False,sep='\t')
del inds, df1_index,df2_index,df1_value,df2_value,ecope,values,i,start,end,G,mir_sim_new,mir_sim_new_matrix
#new-drug
pro_drug_mir_matrix = pro_mir_drug_matrix.T
inds = np.array(pro_drug_mir_matrix.index)
df1_index = np.tile(inds.reshape(-1,1),n).ravel()
df2_index = np.tile(inds.reshape(1,-1),n).ravel()
df1_value = pro_drug_mir_matrix.values.repeat(n,axis=0)
df2_value = np.tile(pro_drug_mir_matrix.values,(n, 1))
ecope = int((n*n)/banch)
values = np.zeros(n*n)
for i in range(ecope):
    start = i * banch
    end = (i+1) * banch
    values[start:end] = (df1_value[start:end] * df2_value[start:end]).sum(axis=1)         
values[ecope * banch: n * n] = (df1_value[ecope * banch: n * n] * df2_value[ecope * banch: n * n]).sum(axis=1) 
G = nx.Graph()
G.add_weighted_edges_from(list(zip(df1_index,df2_index,values)))
drug_sim_new = nx.to_pandas_adjacency(G)
drug_sim_new_matrix = drug_sim_new.loc[drug_name.Drug.tolist(),drug_name.Drug.tolist()]
drug_sim_new_matrix.to_csv('../data/new similarity/Pro-mTD/new-drug-sim.txt',header=None,index=False,sep='\t')
del inds, df1_index,df2_index,df1_value,df2_value,ecope,values,i,start,end,G,drug_sim_new,drug_sim_new_matrix

#读取ncDR_miR-drug关系网络生成关系矩阵
ncDR_assoc_array = pd.read_csv(assoc_path+'ncDR/mir-drug-array.csv')
assoc_data_1 = pd.merge(ncDR_assoc_array,mir_name.miRNA,how='inner',on=['miRNA'])
assoc_data_2 = pd.merge(assoc_data_1,drug_name.Drug,how='inner',on=['Drug'])
assoc_data_3 = pd.merge(assoc_data_2,mir_name.miRNA,how='outer',on=['miRNA'])
assoc_data_3['Drug']=assoc_data_3['Drug'].fillna(value='Cisplatin')
assoc_data_3['Control']=assoc_data_3['Control'].fillna(value=0)
assoc_data_4 = pd.merge(assoc_data_3,drug_name.Drug,how='outer',on=['Drug'])
assoc_data_4['miRNA']=assoc_data_4['miRNA'].fillna(value='hsa-let-7')
assoc_data_4['Control']=assoc_data_4['Control'].fillna(value=0)
del assoc_data_1,assoc_data_2,assoc_data_3,ncDR_assoc_array
temp = list(zip(assoc_data_4.iloc[:,1],assoc_data_4.iloc[:,0],assoc_data_4.iloc[:,2]))
G = nx.Graph()
G.add_weighted_edges_from(temp)
ncDR_mir_drug_matrix = nx.to_pandas_adjacency(G)
ncDR_mir_drug_matrix = ncDR_mir_drug_matrix.loc[mir_name.miRNA.tolist(),drug_name.Drug.tolist()]
ncDR_mir_drug_matrix.to_csv(assoc_path+'ncDR/mir-drug-matrix.csv')
ncDR_mir_drug_matrix.to_csv(assoc_path+'ncDR/mir-drug-matrix.txt',sep='\t',index=None,header=False)
del assoc_data_4,temp,G
#生成新的相似矩阵
m,n = ncDR_mir_drug_matrix.shape
banch = 50000
#new-mirna
inds = np.array(ncDR_mir_drug_matrix.index)
df1_index = np.tile(inds.reshape(-1,1),m).ravel()
df2_index = np.tile(inds.reshape(1,-1),m).ravel()
df1_value = ncDR_mir_drug_matrix.values.repeat(m,axis=0)
df2_value = np.tile(ncDR_mir_drug_matrix.values,(m, 1))
ecope = int((m*m)/banch)
values = np.zeros(m*m)
for i in range(ecope):
    start = i * banch
    end = (i+1) * banch
    values[start:end] = (df1_value[start:end] * df2_value[start:end]).sum(axis=1)         
values[ecope * banch: m * m] = (df1_value[ecope * banch: m * m] * df2_value[ecope * banch: m * m]).sum(axis=1) 
G = nx.Graph()
G.add_weighted_edges_from(list(zip(df1_index,df2_index,values)))
mir_sim_new = nx.to_pandas_adjacency(G)
mir_sim_new_matrix = mir_sim_new.loc[mir_name.miRNA.tolist(),mir_name.miRNA.tolist()]
mir_sim_new_matrix.to_csv('../data/new similarity/ncDR/new-mir-sim.txt',header=None,index=False,sep='\t')
del inds, df1_index,df2_index,df1_value,df2_value,ecope,values,i,start,end,G,mir_sim_new,mir_sim_new_matrix
#new-drug
ncDR_drug_mir_matrix = ncDR_mir_drug_matrix.T
inds = np.array(ncDR_drug_mir_matrix.index)
df1_index = np.tile(inds.reshape(-1,1),n).ravel()
df2_index = np.tile(inds.reshape(1,-1),n).ravel()
df1_value = ncDR_drug_mir_matrix.values.repeat(n,axis=0)
df2_value = np.tile(ncDR_drug_mir_matrix.values,(n, 1))
ecope = int((n*n)/banch)
values = np.zeros(n*n)
for i in range(ecope):
    start = i * banch
    end = (i+1) * banch
    values[start:end] = (df1_value[start:end] * df2_value[start:end]).sum(axis=1)         
values[ecope * banch: n * n] = (df1_value[ecope * banch: n * n] * df2_value[ecope * banch: n * n]).sum(axis=1) 
G = nx.Graph()
G.add_weighted_edges_from(list(zip(df1_index,df2_index,values)))
drug_sim_new = nx.to_pandas_adjacency(G)
drug_sim_new_matrix = drug_sim_new.loc[drug_name.Drug.tolist(),drug_name.Drug.tolist()]
drug_sim_new_matrix.to_csv('../data/new similarity/ncDR/new-drug-sim.txt',header=None,index=False,sep='\t')
del inds, df1_index,df2_index,df1_value,df2_value,ecope,values,i,start,end,G,drug_sim_new,drug_sim_new_matrix




temp = list(zip(pro_assoc_array.iloc[:,0],pro_assoc_array.iloc[:,1],pro_assoc_array.iloc[:,2]))
G = nx.Graph()
G.add_weighted_edges_from(temp)
drug_sim_matrix = nx.to_pandas_adjacency(G)
drug_sim_matrix.to_csv('../data/Drug/drug-sim-matrix.csv')
ncdr_assoc_array = pd.read_csv(assoc_path+'ncDR/mir-drug-array.csv')