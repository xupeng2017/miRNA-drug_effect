#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 09:33:26 2019

@author: spark
"""


import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm

no_edge_name = pd.read_csv('../data/no-edge/ncDR-noedge.csv')
no_edge = no_edge_name.copy()
path = u'../result/KFold-result/ncDR/5-fold/BiRW/array-5fold/'
for i in tqdm(range(1,2501)):
    all_pred = pd.read_csv(path+str(i)+'.csv')
    no_edge_pred = pd.merge(no_edge_name,all_pred,how='left',on=['Drug','miRNA'])
    no_edge['weight-'+str(i)] = no_edge_pred.iloc[:,2]


no_edge['mean'] = no_edge.iloc[:,2:].mean(axis=1)
no_edge_averge = pd.DataFrame({'Drug':no_edge.Drug,'miRNA':no_edge.miRNA,'mean':no_edge.iloc[:,-1]})
no_edge_averge.sort_values(by='mean',ascending=False,inplace=True)
no_edge_averge['Rank'] = range(1,no_edge_averge.shape[0]+1)
no_edge_averge.to_csv('../result/KFold-result/ncDR/5-fold/BiRW/no-edge-rank.csv',index=False)
