# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:23:08 2020

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

#%%

data_2 = pd.read_excel('CopyofJ_LIS.xlsx',skiprows = 2)
Jp = data_2.iloc[:,3]/10000
Jp.index = data_p.columns
Ja = data_2.iloc[:,4]/10000
Ja.index = data_a.columns

#%%
smooth_points = 50
df = Yp
columns_num = (len(df.columns)-1) * smooth_points
index_num =len(df.index)
start_en = df.columns
start_array = np.zeros([index_num , columns_num])
final_df = pd.DataFrame(start_array)
final_df.index = df.index
for j in range(110):
    final_en = np.empty([])
    final_vals = np.empty([])
    
    for i in range(len(df.columns)-1):
        x1 = start_en[i]
        x2 = start_en[i+1]
        y1 = df.iloc[j,i] 
        y2 = df.iloc[j,i+1]
        temp_en = np.linspace(x1,x2,smooth_points)
        b = np.log(y2/y1) / np.log(x2/x1)
        a = y1 / (x1**b)
        
        final_vals_i =  a * temp_en ** b    
        final_vals = np.append(final_vals,final_vals_i)
        final_en = np.append(final_en , temp_en)
        
    final_en = final_en[1:]
    final_vals = final_vals[1:]  
    final_df.iloc[j,:] = final_vals
final_df.columns = final_en
final_df.fillna(0,inplace = True)
Yp = final_df
Yp.fillna(0,inplace = True)

df = Ya
columns_num = (len(df.columns)-1) * smooth_points
index_num =len(df.index)
start_en = df.columns
start_array = np.zeros([index_num , columns_num])
final_df = pd.DataFrame(start_array)
final_df.index = df.index
for j in range(110):
    final_en = np.empty([])
    final_vals = np.empty([])
    
    for i in range(len(df.columns)-1):
        x1 = start_en[i]
        x2 = start_en[i+1]
        y1 = df.iloc[j,i] 
        y2 = df.iloc[j,i+1]
        temp_en = np.linspace(x1,x2,smooth_points)
        b = np.log(y2/y1) / np.log(x2/x1)
        a = y1 / (x1**b)
        
        final_vals_i =  a * temp_en ** b    
        final_vals = np.append(final_vals,final_vals_i)
        final_en = np.append(final_en , temp_en)
        
    final_en = final_en[1:]
    final_vals = final_vals[1:]  
    final_df.iloc[j,:] = final_vals
final_df.columns = final_en
final_df.fillna(0,inplace = True)
Ya = final_df
Ya.fillna(0,inplace = True)
#%%
for i in range(2):
    if i == 0:
        Z = 1
        A = 1
    if i == 1:
        Z = 2
        A = 4
        
    
    Tr = 0.938
    ph = 0.65
    phi = (Z/A) * ph
    
    T = Yp.columns
    Tphi = T + phi
    
    Gamma = (Tr + Tphi) / Tr
    
    Beta = np.sqrt(1 - (1/Gamma)**2)
    
    JLIS = ((2700 * Tphi**1.12) / Beta**2) * ((Tphi +0.67) / 1.67)**-3.93
    if i == 0:
        Jp_cal = JLIS * T * (T + 2*Tr) / Tphi / (Tphi + 2*Tr) /10000
    if i == 1:
        JP_stan = JLIS * T * (T + 2*Tr) / Tphi / (Tphi + 2*Tr) /10000
        Ja_cal = 0.353 * JP_stan #(JLIS * T * (T + 2*Tr) / Tphi / (Tphi + 2*Tr) /10000)

plt.loglog(data_2.iloc[:,0],Jp,'x')
plt.loglog(T,Jp_cal,'b')



#%%

YpxJ = Yp.apply(lambda x: x * Jp_cal,axis = 1)
YpxJ_heights = np.trapz(YpxJ,YpxJ.columns)
Q_p = np.trapz(YpxJ_heights[1:],data_p.index[1:]) + YpxJ_heights[0]

YaxJ = Ya.apply(lambda x: x * Ja_cal,axis = 1)
YaxJ_heights = np.trapz(YaxJ,YaxJ.columns)
Q_a = np.trapz(YaxJ_heights[1:],data_a.index[1:]) + YaxJ_heights[0]
