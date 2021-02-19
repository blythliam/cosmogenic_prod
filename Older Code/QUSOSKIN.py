# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 12:52:01 2020

@author: lpb20
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import pi





#%%

data_p = pd.read_excel('COproduction.xlsx', '14C_p',skiprows=4,index_col=0)
data_p.fillna(0,inplace = True)
data_a = pd.read_excel('COproduction.xlsx', '14C_a',skiprows=4,index_col=0)
data_a.fillna(0,inplace = True)


Yp = data_p * pi
Ya = data_a * pi

data_2 = pd.read_excel('CopyofJ_LIS.xlsx',skiprows = 2)
Jp = data_2.iloc[:,3]/10000
Jp.index = data_p.columns
Ja = data_2.iloc[:,4]/10000
Ja.index = data_a.columns






YpxJ = Yp.apply(lambda x: x * Jp,axis = 1)
YpxJ_heights = np.trapz(YpxJ,YpxJ.columns)
Q_p = np.trapz(YpxJ_heights[1:],data_p.index[1:]) + YpxJ_heights[0]

YaxJ = Ya.apply(lambda x: x * Ja,axis = 1)
YaxJ_heights = np.trapz(YaxJ,YaxJ.columns)
Q_a = np.trapz(YaxJ_heights[1:],data_a.index[1:]) + YaxJ_heights[0]
 

#%%













