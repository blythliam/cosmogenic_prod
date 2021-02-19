# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 14:13:52 2020

@author: lpb20
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import pi


#%% Polar production


data_p = pd.read_excel('COproduction.xlsx', '14C_p',skiprows=4,index_col=0)
data_p.fillna(0,inplace = True)
data_a = pd.read_excel('COproduction.xlsx', '14C_a',skiprows=4,index_col=0)
data_a.fillna(0,inplace = True)
M = 8
phi = 0.65
T_r = 0.938
#%% Smoothing data
smooth_points = 50

df = data_p
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

data_p = final_df

df = data_a
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

data_a = final_df

#plt.plot(final_df.columns, final_df.iloc[i,:],'r')
#plt.plot(data_p.columns, data_p.iloc[i,:],'b')





#%% Producing Q for protons

Z = 1
A = 1
phi_o = phi * (Z/A)
T = np.array(data_p.columns) #Using the energies from data not own energies
#    T = np.linspace(0,1000,10000)
Tphi = T + phi_o
PTphi = (Tphi*(Tphi+2*T_r))**0.5
JLIS = (1.9e4 * PTphi**-2.78) / (1 + 0.4866 * PTphi**-2.51)
JP = JLIS * T * (T + 2 * T_r) / (T + phi_o) / (T + phi_o + 2*T_r)/10000

Yp = data_p * pi    
YpxJ = Yp * JP
YpxJ_int = YpxJ.apply(lambda x: np.trapz(x,YpxJ.columns),axis = 1)

df = YpxJ
I_height = []
for j in range(110):
    df = YpxJ.iloc[j,:]
    I = 0

    for i in range(22):
        x1 = df.index[i]    
        x2 = df.index[i+1]
        y1 = df.values[i]
        y2 = df.values[i+1]
        if y1 == 0 or y2 ==0:
            I_now = 0
        else:
            print(j)
            b = np.log(y2/y1) / np.log(x2/x1)
            a = y1 / (x1 **b)
            I_now = (a/(b+1)) * (x2**(b+1) - x1**(b+1)) #* (x2 - x1)
     #   I_now = (x2 - x1)*y1  +  (a/(b+1))*(x2**(b+1) * x1**(b+1))
        I+=I_now
    I_height.append(I)


 
Qp = np.trapz(YpxJ_int[1:], YpxJ.index[1:]) + YpxJ_int.iloc[0]


#%% Producing J for alpha
Z = 2
A = 4
phi_o = phi * (Z/A)
T = np.array(data_a.columns) #Using the energies from data not own energies
#    T = np.linspace(0,1000,10000)
Tphi = T + phi_o
PTphi = (Tphi*(Tphi+2*T_r))**0.5

JLIS = 0.05 * (1.9e4 * PTphi**-2.78) / (1 + 0.4866 * PTphi**-2.51)  
JA = JLIS * T * (T + 2 * T_r) / (T + phi_o) / (T + phi_o + 2*T_r)/10000

Ya = data_a * pi
YaxJ = Ya * JA
YaxJ.columns = YaxJ.columns * 4
YaxJ_int = YaxJ.apply(lambda x: np.trapz(x,YaxJ.columns),axis = 1)
Qa = np.trapz(YaxJ_int[1:], YaxJ.index[1:]) + YaxJ_int.iloc[0]

Q_tot = Qp + Qa
    
#%%
# Smothing curve

original = Yp.iloc[0,:].values
energies_orig = Yp.iloc[0].index

new_orig= []
new_energies= []
for i in range(22):
    x1 = energies_orig[i]
    x2 = energies_orig[i+1]
    y1 = original[i]
    y2 = original[i+1]
    new_ens = np.linspace(x1,x2,1000)
    b = np.log(y2/y1) / np.log(x2/x1)
    a = y1 / (x1**b)
    new_ys = a * new_ens**b
    new_orig.append(new_ys)
    new_energies.append(new_ens)
new_vals = np.array(new_orig).reshape([22000,1])
new_energies = np.array(new_energies).reshape([22000,1])
#plt.plot(new_energies2,new_vals)

Z = 1
A = 1
phi_o = phi * (Z/A)
T = np.array(new_energies) #Using the energies from data not own energies
#    T = np.linspace(0,1000,10000)
Tphi = T + phi_o
PTphi = (Tphi*(Tphi+2*T_r))**0.5
JLIS = (1.9e4 * PTphi**-2.78) / (1 + 0.4866 * PTphi**-2.51)
JP = JLIS * T * (T + 2 * T_r) / (T + phi_o) / (T + phi_o + 2*T_r)/10000

YxJ = JP * new_vals

plt.plot(new_energies, YxJ)
plt.plot(YpxJ.iloc[0,:].index, YpxJ.iloc[0,:].values)


#%% Creating a smoothed dataframe

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


#%% Using proper integral







df = Yp.iloc[0,:]

I = 0

for i in range(22):
    x1 = df.index[i]    
    x2 = df.index[i+1]
    y1 = df.values[i]
    y2 = df.values[i+1]
    b = np.log(y2/y1) / np.log(x2/x1)
    a = y1 / (x1 **b)
    I_now = (a/(b+1)) * (x2**(b+1) - x1**(b+1)) #* (x2 - x1)
 #   I_now = (x2 - x1)*y1  +  (a/(b+1))*(x2**(b+1) * x1**(b+1))
    print(I_now)
    I+=I_now




















    
    
    
    