# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 23:28:00 2019

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

drugname = pd.read_csv('../data/Drug/ncDR-database/DrugName.txt',header=None)
mirname = pd.read_csv('../data/miRNA/miRName.txt',header=None)
mir_drug_ass_array = pd.read_csv('../data/ncDR-database/mir-drug-array.csv')
path = u'../result/5-fold/ncDR-database/BiRW/'
path_kf = u'../data/Kfold/5-fold/ncDR-database/'
result_5fold = []

for i in tqdm(range(1, 501)):
    #读取mir-drug-pred矩阵
    drug_mir_pred_matrix = pd.read_csv(path+'matrix-5fold/'+str(i)+'.csv',header=None)
#    drug_mir_pred_matrix.index = drugname.iloc[:,0]
#    drug_mir_pred_matrix.columns = mirname.iloc[:,0]
    m,n = drug_mir_pred_matrix.shape
    inds = np.array(drugname.iloc[:,0])
    cols = np.array(mirname.iloc[:,0]) 
    Drug = np.tile(inds.reshape(-1, 1), n).ravel()
    miRNA = np.tile(cols.reshape(1, -1), m).ravel()
    weight = drug_mir_pred_matrix.values.ravel()
    drug_mir_array = pd.DataFrame({'Drug':Drug,'miRNA':miRNA,'Weight':weight})
    drug_mir_pred_array = drug_mir_array.sort_values(by='Weight',ascending=False)
    drug_mir_pred_array.to_csv(path+'array-5fold/'+str(i)+'.csv',index=None)
    del drug_mir_pred_matrix,m,n,inds,cols,Drug,miRNA,weight,drug_mir_array
    
    test_data = pd.read_csv(path_kf+'test/'+str(i)+'.csv')
#    train_data = pd.read_csv(path_kf+'train/'+str(i)+'.csv')
    #筛选有边无边的预测值
    test_pred = pd.merge(test_data.iloc[:,0:2],drug_mir_pred_array,how='inner',on=['Drug','miRNA'])
#    train_pred = pd.merge(train_data.iloc[:,0:2],drug_mir_pred_array,how='inner',on=['Drug','miRNA'])
    dele_pred = pd.merge(mir_drug_ass_array.iloc[:,0:2],drug_mir_pred_array,how='inner',on=['Drug','miRNA'])
    df1 = pd.concat([drug_mir_pred_array,dele_pred],ignore_index=True,sort=False)
    noedge_pred = df1.drop_duplicates(keep=False)
    #添加正负标签
    test_pred['label']=1
    col_name=noedge_pred.columns.tolist() 
    col_name.insert(3,'label')
    df3=noedge_pred.reindex(columns=col_name)
    df3['label'] = 0
    #df3['label']=df3['label'].fillna(value=0)
    #拼接正负样本
    pred = pd.concat([test_pred,df3],ignore_index=True)
    y_true = pred['label'].values
    y_scores = pred['Weight'].values
    #计算AUC
    AUC = roc_auc_score(y_true, y_scores)
    result_5fold.append(AUC)

#保存AUC值
result = pd.DataFrame({'5fold-auc':result_5fold})
result.to_csv('../result/5-fold/ncDR-database/BiRW/5fold-auc.csv',index=None)
    
    
