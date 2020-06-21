# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 15:29:31 2019

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

pro_mtd_array = pd.read_csv('../data/miRNA-Drug association/Pro-mTD/mir-drug-association.csv')
ncdr_array = pd.read_csv('../data/miRNA-Drug association/ncDR/mir-drug-array.csv')
#两个数据集公共关系
merge_data = pd.merge(pro_mtd_array,ncdr_array,how='inner',on=['miRNA','Drug'])
merge_data.iloc[:,0:2].to_csv('../result/merge-array.csv',index=None)
#预测结果与原始关系的公共部分
pro_pred = pd.read_csv('../result/KFold-result/Pro-mTD/5-fold/BiRW/no-edge-rank.csv')
ncdr_pred = pd.read_csv('../result/KFold-result/ncDR/5-fold/BiRW/no-edge-rank.csv')
pro_rank = pro_pred.iloc[0:500,:]
ncdr_rank = ncdr_pred.iloc[0:500,:]
#保存预测排名的top-500
pro_rank.to_csv('../result/analysis/pro-mTD-predict-top500.csv',index=None)
ncdr_rank.to_csv('../result/analysis/ncDR-predict-top500.csv',index=None)
#找预测结果与另外关系网络的公共部分
pro_pred_ncdr_merge = pd.merge(pro_rank,ncdr_array.iloc[:,0:2],how='inner',on=['miRNA','Drug'])
pro_pred_ncdr_merge.to_csv('../result/analysis/pro-mTD-top500-ncDR-array-share.csv',index=None)
ncdr_pred_pro_merge = pd.merge(ncdr_rank,pro_mtd_array.iloc[:,0:2],how='inner',on=['miRNA','Drug'])
ncdr_pred_pro_merge.to_csv('../result/analysis/ncDR-top500-pro-mTD-array-share.csv',index=None)
pro_ncdr_pred_share = pd.merge(pro_rank,ncdr_rank,how='inner',on=['miRNA','Drug'])
pro_ncdr_pred_share.to_csv('../result/analysis/pro-ncDR-top500-share-array.csv',index=None)
