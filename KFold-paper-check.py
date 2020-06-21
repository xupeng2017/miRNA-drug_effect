# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 19:35:42 2019

@author: Hp
"""

import numpy as np
import pandas as pd
#读取预测结果
pro_rank = pd.read_excel('../result/analysis/top-rank-paper.xlsx',sheet_name='Pro-mTD-paper')
ncdr_rank = pd.read_excel('../result/analysis/top-rank-paper.xlsx',sheet_name='ncDR-paper')
df = pd.read_excel('../../mir-drug-sample/result/ori-result/no-edge-paper.xlsx',sheet_name='Sheet1')
df1 = pd.DataFrame({'Drug':df.Drug,'miRNA':df.miRNA,'PMID':df.PMID,'mir-col':df['mir-col'],'drug-col':df['drug-col']})
df2 = pd.DataFrame({'Drug':pro_rank.Drug,'miRNA':pro_rank.miRNA,'PMID':pro_rank.PMID,'mir-col':pro_rank['mir-col'],'drug-col':pro_rank['drug-col']})
df3 = pd.merge(pro_rank,df1,how='left',on=['Drug','miRNA'])
df3.to_csv('../result/df3.csv')
df4 = pd.read_csv('../result/df3.csv')
df5 = pd.merge(df4,df2,how='left',on=['Drug','miRNA'])
