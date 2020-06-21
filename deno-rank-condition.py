# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 20:18:19 2019

@author: Hp
"""

import numpy as np
import pandas as pd
from tqdm import tqdm


df = pd.DataFrame()
#读取rank结果
for i in range(1,99):
    rank_result = pd.read_csv('../result/deno/validate-rank/Pro-mTD/BiRW/'+str(i)+'.csv')
    df = pd.concat([df,rank_result],ignore_index= False)
df.to_csv('../result/deno/validate-rank/Pro-mTD/Pro-mTD-validate-BiRW-rank.csv',index=None)
