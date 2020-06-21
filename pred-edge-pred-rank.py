# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 09:41:49 2019

@author: Hp
"""

import os 
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm

#miRNA-Drug关系数组
pro_assoc_array = pd.read_csv('../data/miRNA-Drug association/Pro-mTD/mir-drug-association.csv')
ncdr_assoc_array = pd.read_csv('../data/miRNA-Drug association/ncDR/mir-drug-array.csv')
#miRNA-Drug预测数组
pro_bi_pred_array = pd.read_csv('../result/Pro-mTD/BiRW/mir-drug-predict-array.csv',index_col=None)
ncdr_bi_pred_array = pd.read_csv('../result/ncDR/BiRW/mir-drug-predict-array.csv',index_col=None)
pro_bi_pred_array.sort_values(by='Weight',ascending=False,inplace=True)
ncdr_bi_pred_array.sort_values(by='Weight',ascending=False,inplace=True)
pro_bi_pred_array['Rank'] = range(1,pro_bi_pred_array.shape[0]+1)
ncdr_bi_pred_array['Rank'] = range(1,ncdr_bi_pred_array.shape[0]+1)
#无边值的预测分布
pro_bi_edge = pd.merge(pro_bi_pred_array,pro_assoc_array.iloc[:,0:2],how='inner',on=['miRNA','Drug'])
pro_bi_no = (pro_bi_pred_array.append(pro_bi_edge,sort=False)).drop_duplicates(keep=False)
pro_bi_no['Rank'] = pro_bi_no.Rank-pro_assoc_array.shape[0]
pro_bi_no.to_csv('../result/Pro-mTD/BiRW/no-edges-predict.csv',index=None)
ncdr_bi_edge = pd.merge(ncdr_bi_pred_array,ncdr_assoc_array.iloc[:,0:2],how='inner',on=['miRNA','Drug'])
ncdr_bi_no = (ncdr_bi_pred_array.append(ncdr_bi_edge,sort=False)).drop_duplicates(keep=False)
ncdr_bi_no['Rank'] = ncdr_bi_no.Rank-ncdr_assoc_array.shape[0]
ncdr_bi_no.to_csv('../result/ncDR/BiRW/no-edges-predict.csv',index=None)
