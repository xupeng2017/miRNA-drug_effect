# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 09:27:19 2019

@author: Hp
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import networkx as nx

pro_birw = pd.read_csv('../result/KFold-result/Pro-mTD/5-fold/BiRW/5fold-auc.csv')
pro_rw = pd.read_csv('../result/KFold-result/Pro-mTD/5-fold/RW/5fold-auc.csv')
ncdr_birw = pd.read_csv('../result/KFold-result/ncDR/5-fold/BiRW/5fold-auc.csv')
ncdr_rw = pd.read_csv('../result/KFold-result/ncDR/5-fold/RW/5fold-auc.csv')
pro_birw['Database'] = 'Pro-mTD'
pro_birw['Method'] = 'BiRW'
pro_rw['Database'] = 'Pro-mTD'
pro_rw['Method'] = 'RW'
ncdr_birw['Database'] = 'ncDR'
ncdr_birw['Method'] = 'BiRW'
ncdr_rw['Database'] = 'ncDR'
ncdr_rw['Method'] = 'RW'
#auc拼接
auc_concat = pd.concat([pro_birw,pro_rw,ncdr_birw,ncdr_rw])
pal_style = ['deep', 'muted', 'pastel', 'bright', 'dark','colorblind']

plt.figure(figsize=(10,8))
ax = sns.boxplot(x="Database", y="AUC",hue="Method",data=auc_concat,palette="bright")
plt.ylim(0.95,1.0)
sns.palplot(sns.color_palette("hls", 8))
