# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 22:40:13 2019

@author: Hp
"""
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import networkx as nx

rank_result = pd.read_excel('../result/analysis/deno-top-statistic.xlsx',sheet_name='Sheet1')
pro_rank = rank_result.iloc[0:10,0:4]
ncdr_rank = rank_result.iloc[10:20,0:4]
#画图
plt.style.use('ggplot')
plt.figure(12,figsize=(14,8))
ax1 = plt.subplot(121)
sns.barplot(x='Top',y='Count',hue='Method',data=pro_rank)
ax1.set_title('Pro-mTD Dataset')
ax2 = plt.subplot(122)
sns.barplot(x='Top',y='Count',hue='Method',data=ncdr_rank)
ax2.set_title('ncDR Dataset')
plt.tight_layout()
plt.show()