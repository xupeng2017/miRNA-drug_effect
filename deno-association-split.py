# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 19:22:11 2019

@author: Hp
"""

import os 
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import auc, roc_curve
from sklearn.metrics import roc_auc_score
from sklearn import metrics
#原始关系数组
pro_mir_drug = pd.read_csv('../data/miRNA-Drug association/Pro-mTD/mir-drug-association.csv')
#关系网络药物个数
drug_list = pro_mir_drug[['Drug']].drop_duplicates(keep='first')
banch = 50000
#读取药物和miRNA名
drug_name = pd.read_csv('../data/Drug/DrugName.txt',sep='\t',header=None)
new_col = ['Drug']
drug_name.columns= new_col
mir_name = pd.read_csv('../data/miRNA/miRName.txt',sep='\t',header=None)
new_col = ['miRNA']
mir_name.columns= new_col
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

    
for i in range(drug_list.shape[0]):
    very = pro_mir_drug[pro_mir_drug.Drug!=drug_list.iloc[i,0]]
    vary_rank_arry = pro_mir_drug[pro_mir_drug.Drug==drug_list.iloc[i,0]]
    vary_rank_arry.to_csv('../data/deno/Pro-mTD/new-association/'+str(i)+'.csv',index=None)
    #生成关系矩阵
    df1 = pd.merge(very,mir_name.miRNA,how='inner',on=['miRNA'])
    df2 = pd.merge(df1,drug_name.Drug,how='inner',on=['Drug'])
    df3 = pd.merge(df2,mir_name.miRNA,how='outer',on=['miRNA'])
    df3['Drug']=df3['Drug'].fillna(value='Doxorubicin')
    df3['Control']=df3['Control'].fillna(value=0)
    df4 = pd.merge(df3,drug_name.Drug,how='outer',on=['Drug'])
    df4['miRNA']=df4['miRNA'].fillna(value='hsa-let-7')
    df4['Control']=df4['Control'].fillna(value=0)
    mir_drug_matrix = array_matrix(df4)
    mir_drug_matrix.to_csv('../data/deno/Pro-mTD/matrix/'+str(i)+'.txt',sep='\t',index=None,header=False)
    #生成新的相似矩阵
    #生成新的miRNA相似矩阵
    new_mir_sim = new_simliar(mir_drug_matrix)
    new_mir_sim_matrix = new_mir_sim.loc[mir_name.miRNA.tolist(),mir_name.miRNA.tolist()]
    new_mir_sim_matrix.to_csv('../data/deno/Pro-mTD/new-mir/'+str(i)+'.txt',sep='\t',index=None,header=False)
    #生成新的drug相似矩阵
    new_drug_sim = new_simliar(mir_drug_matrix.T)
    new_drug_sim_matrix = new_drug_sim.loc[drug_name.Drug.tolist(),drug_name.Drug.tolist()]
    new_drug_sim_matrix.to_csv('../data/deno/Pro-mTD/new-drug/'+str(i)+'.txt',sep='\t',index=None,header=False)
    
    
    
