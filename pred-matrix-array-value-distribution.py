# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 22:31:55 2019

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

#drug、miRNA名
drug_name = pd.read_csv('../data/Drug/DrugName.txt',sep='\t',header=None)
mir_name = pd.read_csv('../data/miRNA/miRName.txt',sep='\t',header=None)
#miRNA-Drug关系数组
pro_assoc_array = pd.read_csv('../data/miRNA-Drug association/Pro-mTD/mir-drug-association.csv')
ncdr_assoc_array = pd.read_csv('../data/miRNA-Drug association/ncDR/mir-drug-array.csv')
#Bi-RW和RW预测矩阵
pro_bi_pred_matrix = pd.read_csv('../result/Pro-mTD/BiRW/mir-drug-predict-matrix.csv',header=None)
pro_rw_pred_matrix = pd.read_csv('../result/Pro-mTD/RW/mir-drug-predict-matrix.csv',header=None)
ncdr_bi_pred_matrix = pd.read_csv('../result/ncDR/BiRW/mir-drug-predict-matrix.csv',header=None)
ncdr_rw_pred_matrix = pd.read_csv('../result/ncDR/RW/mir-drug-predict-matrix.csv',header=None)

#转为数组
inds = np.array(drug_name.iloc[:,0])
cols = np.array(mir_name.iloc[:,0]) 

def matrix_array(df):
    df = pd.DataFrame(df)
    m,n = df.shape
    Drug = np.tile(inds.reshape(-1, 1), n).ravel()
    miRNA = np.tile(cols.reshape(1, -1), m).ravel()
    weight = df.values.ravel()
    result_array = pd.DataFrame({'Drug':Drug,'miRNA':miRNA,'Weight':weight})
    result_array.sort_values(by='Weight',ascending=False,inplace=True)
    return result_array
#保存预测结果数组
pro_bi_pred_array = matrix_array(pro_bi_pred_matrix)
pro_bi_pred_array.to_csv('../result/Pro-mTD/BiRW/mir-drug-predict-array.csv',index=None)

pro_rw_pred_array = matrix_array(pro_rw_pred_matrix)
pro_rw_pred_array.to_csv('../result/Pro-mTD/RW/mir-drug-predict-array.csv',index=None)

ncdr_bi_pred_array = matrix_array(ncdr_bi_pred_matrix)
ncdr_bi_pred_array.to_csv('../result/ncDR/BiRW/mir-drug-predict-array.csv',index=None)

ncdr_rw_pred_array = matrix_array(ncdr_rw_pred_matrix)
ncdr_rw_pred_array.to_csv('../result/ncDR/RW/mir-drug-predict-array.csv',index=None)

#有边、无边预测值分布
#pro-birw-predict
have_pro_bi = pd.merge(pro_bi_pred_array,pro_assoc_array.iloc[:,0:2],how='inner',on=['Drug','miRNA'])
no_pro_bi = (pro_bi_pred_array.append(have_pro_bi,sort=False)).drop_duplicates(keep=False)
y1 = pd.DataFrame({'Weight':have_pro_bi.Weight,'Edges':'yes','database-method':'Pro-mTD_BiRW'})
y2 = pd.DataFrame({'Weight':no_pro_bi.Weight,'Edges':'no','database-method':'Pro-mTD_BiRW'})
#pro-rw-predict
have_pro_rw = pd.merge(pro_rw_pred_array,pro_assoc_array.iloc[:,0:2],how='inner',on=['Drug','miRNA'])
no_pro_rw = (pro_rw_pred_array.append(have_pro_rw,sort=False)).drop_duplicates(keep=False)
y3 = pd.DataFrame({'Weight':have_pro_rw.Weight,'Edges':'yes','database-method':'Pro-mTD_RW'})
y4 = pd.DataFrame({'Weight':no_pro_rw.Weight,'Edges':'no','database-method':'Pro-mTD_RW'})
#ncDR-birw-predict
have_ncdr_bi = pd.merge(ncdr_bi_pred_array,ncdr_assoc_array.iloc[:,0:2],how='inner',on=['Drug','miRNA'])
no_ncdr_bi = (ncdr_bi_pred_array.append(have_ncdr_bi,sort=False)).drop_duplicates(keep=False)
y5 = pd.DataFrame({'Weight':have_ncdr_bi.Weight,'Edges':'yes','database-method':'ncDR_BiRW'})
y6 = pd.DataFrame({'Weight':no_ncdr_bi.Weight,'Edges':'no','database-method':'ncDR_BiRW'})
#ncDR-rw-predict
have_ncdr_rw = pd.merge(ncdr_rw_pred_array,ncdr_assoc_array.iloc[:,0:2],how='inner',on=['Drug','miRNA'])
no_ncdr_rw = (ncdr_rw_pred_array.append(have_ncdr_rw,sort=False)).drop_duplicates(keep=False)
y7 = pd.DataFrame({'Weight':have_ncdr_rw.Weight,'Edges':'yes','database-method':'ncDR_RW'})
y8 = pd.DataFrame({'Weight':no_ncdr_rw.Weight,'Edges':'no','database-method':'ncDR_RW'})

pred_concat = pd.concat([y1,y2,y3,y4,y5,y6,y7,y8,])
plt.figure(figsize=(16,10))
ax = sns.boxplot(x="database-method", y="Weight", hue="Edges",data=pred_concat, palette="Set3")
plt.show()





