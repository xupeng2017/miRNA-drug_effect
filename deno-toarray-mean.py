# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 10:09:30 2019

@author: Hp
"""

import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import auc, roc_curve
from sklearn.metrics import roc_auc_score
from sklearn import metrics
import os

#miRNA和Drug名
drug_name = pd.read_csv('../data/Drug/DrugName.txt',sep='\t',header=None)
mir_name = pd.read_csv('../data/miRNA/miRName.txt',sep='\t',header=None)

#将矩阵转为数组
def matrix_array(df):
    df = pd.DataFrame(df)
    m,n = df.shape
    inds = np.array(drug_name.iloc[:,0])
    cols = np.array(mir_name.iloc[:,0]) 
    Drug = np.tile(inds.reshape(-1, 1), n).ravel()
    miRNA = np.tile(cols.reshape(1, -1), m).ravel()
    weight = df.values.ravel()
    df_array = pd.DataFrame({'Drug':Drug,'miRNA':miRNA,'Weight':weight})
    return df_array


for i in tqdm(range(1,69)):
    validate_array = pd.read_csv('../data/deno/ncDR/new-association/'+str(i)+'.csv')
    a = validate_array.iloc[0,1]
    #读取预测矩阵
    pred_matrix = pd.read_csv('../result/deno/ncDR/RW/matrix/'+str(i)+'.csv',header=None)
    pred_array = matrix_array(pred_matrix)
    #筛选需要的边
    no_edge_pred = pred_array[pred_array.Drug==a]
    df = no_edge_pred.copy()
    df.sort_values(by='Weight',ascending=False,inplace=True)
    df['Rank'] = range(1,df.shape[0]+1)
    validate_rank = pd.merge(validate_array.iloc[:,0:2],df,how='inner',on=['Drug','miRNA'])
    validate_rank.sort_values(by='Rank',ascending=True,inplace=True)
    validate_rank.to_csv('../result/deno/validate-rank/ncDR/RW/'+str(i)+'.csv',index=None)
    
  