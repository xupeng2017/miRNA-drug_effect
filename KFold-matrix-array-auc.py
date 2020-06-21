# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 22:43:12 2019

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

#miRNA、Drug-Name
drug_name = pd.read_csv('../data/Drug/DrugName.txt',sep='\t',header=None)
mir_name = pd.read_csv('../data/miRNA/miRName.txt',sep='\t',header=None)
#mir-drug assocation
mir_drug_array = pd.read_csv('../data/miRNA-Drug association/ncDR/mir-drug-array.csv')
#no-edge pair name
no_edge_pair = pd.read_csv('../data/no-edge/ncDR-noedge.csv')

train_path = u'../data/KFold/5-fold/ncDR/train/'
pred_path = u'../result/KFold-result/ncDR/5-fold/'

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
    
result_5fold = []


for i in tqdm(range(1, 2501)):
    #train集的关系
    train_array = pd.read_csv(train_path+str(i)+'.csv')
    #Drig-miRNA预测结果矩阵
    drug_mir_pred_matrix = pd.read_csv(pred_path+'BiRW/matrix-5fold/'+str(i)+'.csv',header=None)
    #转为数组
    pred_array = matrix_array(drug_mir_pred_matrix)
    pred_array.to_csv(pred_path+'BiRW/array-5fold/'+str(i)+'.csv',index=None)
    #筛选有边无边预测结果
    dele_pred = pd.merge(pred_array,mir_drug_array.iloc[:,0:2],how='inner',on=['Drug','miRNA'])
    noedge_pred = (pred_array.append(dele_pred,sort=False)).drop_duplicates(keep=False)
    train_pred = pd.merge(dele_pred,train_array.iloc[:,0:2],how='inner',on=['Drug','miRNA'])
    havedge_pred = (dele_pred.append(train_pred,sort=False)).drop_duplicates(keep=False)
    #保存无边的预测值
    #no_edge_pred_list = pd.merge(no_edge_pair,noedge_pred,how='inner',on=['Drug','miRNA'])
    #no_edge_pair['weight-'+str(i)] =no_edge_pred_list.iloc[:,i+1]
    
    
    #添加标签
    col_name=havedge_pred.columns.tolist() 
    col_name.insert(3,'label')
    true_sample = havedge_pred.reindex(columns=col_name)
    true_sample['label'] = 1
    false_sample = noedge_pred.reindex(columns=col_name)
    false_sample['label'] = 0
    #拼接预测结果
    pred_result = pd.concat([true_sample,false_sample],ignore_index=True)
    y_true = pred_result['label'].values
    y_scores = pred_result['Weight'].values
    #计算AUC
    AUC = roc_auc_score(y_true, y_scores)
    result_5fold.append(AUC)
    del train_array,drug_mir_pred_matrix,pred_array,dele_pred,noedge_pred,train_pred,havedge_pred,col_name,true_sample,false_sample,pred_result,AUC


result = pd.DataFrame({'AUC':result_5fold,'Database-Method':'ncDR_BiRW'})
result.to_csv(pred_path+'BiRW/5fold-auc.csv',index=None)
#df = no_edge_pred_list.iloc[:,2:].mean(axis=1)
#no_edge_pair['no_edge_mean'] = df.apply(lambda x:x.sum()/2500,axis=1)