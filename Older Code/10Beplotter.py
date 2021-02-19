# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 10:31:24 2020

@author: lpb20
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 10:24:05 2020

@author: lpb20

An attempt to recreate figure 8
"""

def smooth_df(df,smooth_points = 50):
    from scipy import interpolate
    import pandas as pd
    energies = np.log(df.index)
#    print(len(energies))
    new_energies = np.linspace(energies[0],energies[-1],len(energies)*10)
    vals = np.log(df.fillna(0).values)
#    print(len(vals))
    f = interpolate.interp1d(energies, vals,kind = 'linear')
    new_values = np.exp(f(new_energies))
    new_values = pd.Series(new_values)
    new_values.index = np.exp(new_energies)
    return new_values



def LIS(E_vals, phi, p_or_a):
    Tr = 0.938
    if p_or_a == 'p':
        Z = 1
        A = 1
        phi = (Z/A) * phi
        T = E_vals
        Tphi = T + phi
        Gamma = (Tr + Tphi) / Tr
        Beta = np.sqrt(1 - (1/Gamma)**2)
        JLIS = ((0.27 * Tphi**1.12) / Beta**2) * ((Tphi +0.67) / 1.67)**-3.93
        J = JLIS * T * (T + 2*Tr) / Tphi / (Tphi + 2*Tr)
        J = pd.Series(J)
        J.index = E_vals
        return J
    
    if p_or_a == 'a':
        Z = 2
        A = 4
        phi = (Z/A) * phi
        T = E_vals
        Tphi = T + phi
        Gamma = (Tr + Tphi) / Tr
        Beta = np.sqrt(1 - (1/Gamma)**2)
        JLIS = ((0.27 * Tphi**1.12) / Beta**2) * ((Tphi +0.67) / 1.67)**-3.93
        J = 0.3 * JLIS * T * (T + 2*Tr) / Tphi / (Tphi + 2*Tr)
        J.index = E_vals
        return J

def Gen_PC():  
    R0 = 6371.2e3
    
    g10 =-29992e-9
    g11 = -1956e-9
    h11 = 5604e-9
    g20 = -1997e-9
    g21 = 3027e-9
    h21 = -2129e-9
    h22 = -200e-9
    g22 = 1663e-9
    
    B0 = np.sqrt(g10**2 + g11**2 + h11**2)
    M = (4*pi / 4*pi*1e-7) * B0 * R0 **3
    
    L0 = 2* g10 * g20 + np.sqrt(3) * (g11*g21 + h11*h21)
    L1 = -1 * g10 * g20 + np.sqrt(3) * (g10*g21 + g11*g22 + h11*h22)
    L2 = -1 * h11 * g20 + np.sqrt(3) * (g10*h21 - h11*g22 + g11*h22)
    
    E = (L0*g10 + L1 * g11  + L2 * h11) / 4 * B0**2
    
    x0 = ((L1-g11*E) / 3*B0**2) * R0
    y0 = ((L2-h11*E) / 3*B0**2) * R0
    z0 = ((L0-g10*E) / 3*B0**2) * R0
    
    
    
    lats = np.linspace(0, pi, 100)
    longs = np.linspace(-pi, pi, 100)
    
    long, lat = np.meshgrid(longs,lats)
    d = np.sqrt(x0**2 + y0**2 + z0**2)
    R2 = R0**2 + d - 2 * R0 * (x0 * np.sin(lat) * np.cos(long) + y0 * np.sin(lat) * np.sin(long) + z0*np.cos(lat)) 
    R = np.sqrt(np.abs(R2))
    cosphiG = (1/(R * B0)) * (g11*(x0 - R0 * np.sin(lat) * np.cos(long)) + h11 * (y0 - R0 * np.sin(lat) * np.sin(long)) + g10 * (z0 - R0 * np.cos(lat)))
    
    phiG = np.arccos(cosphiG)
    
    Pc = 1.9 * M * (R0 / R)**2 * np.sin(phiG)**4 / 1e9
    
    return Pc
    
    
    
    
    
    




import pandas as pd
import cairo
import matplotlib.pyplot as plt
#%matplotlib qt
import numpy as np
from math import pi

#Deciding whether to print plots or not
plot = False
M= 7.8
Er = 0.938
#%%Reading in data


S_p = pd.read_excel('COproduction.xlsx', sheet_name='14C_p', skiprows=4, index_col=0)
S_a = pd.read_excel('COproduction.xlsx', sheet_name='14C_a', skiprows=4, index_col=0)


phi_data = pd.read_csv('phi_over_time.txt',sep = '  ',skiprows = 9, header=None,engine=('python'))


Sp_smooth = pd.DataFrame(np.zeros([110,230]))
Sa_smooth = pd.DataFrame(np.zeros([110,230]))


for i in range(110):
    Sp_0_smooth = smooth_df(S_p.iloc[i])
    Sp_smooth.iloc[i] = Sp_0_smooth.values
    Sa_0_smooth = smooth_df(S_a.iloc[i])
    Sa_smooth.iloc[i] = Sa_0_smooth.values
    
Sp_smooth.index = S_p.index
Sp_smooth.columns = Sp_0_smooth.index
Sa_smooth.index = S_a.index
Sa_smooth.columns = Sa_0_smooth.index
#Plotting Source function as function of altidude for specific energy
#data = pd.read_csv('Clproduction.csv', skiprows=4, index_col=0).iloc[:-1]
#    S.columns = S.columns.astype(float)
#    S.fillna(0,inplace = True)
Y_p = pi * Sp_smooth # Y is the smoothed version of S * pi 
Y_a = pi * Sa_smooth

phis = np.linspace(0,2,20)
phis_glob = []
for phi in phis:
    print(phi)
    J_p = LIS(Y_p.columns, phi, 'p') 
    J_a = LIS(Y_a.columns, phi, 'a')
    
    YpxJp = Y_p.apply(lambda x: x*J_p,axis = 1)
    YpxJa = Y_a.apply(lambda x: x*J_a,axis = 1)
    
    
    
    lats = np.linspace( -pi/2 , pi/2 , 180)
    lat_int = []
    for lat in lats:
#        print(lat)
        Pc = 1.9 * M * (np.cos(lat))**4
        Ecp = 0.938 * (np.sqrt(1+((Pc)/(0.938))**2)-1)
        Eca = 0.938 * (np.sqrt(1 + ((2 * Pc)/(4 *0.938))**2) -1)
        
        YpxJa_above = YpxJa[YpxJa.columns[YpxJa.columns>Eca]]
        
        YpxJp_above = YpxJp[YpxJp.columns[YpxJp.columns>Ecp]]
        
        
        def Int(y,x):
            I = np.sum(np.nan_to_num(np.diff(y*x) / ((np.diff(np.log(y))/np.diff(np.log(x))) + 1),nan=0.0))
            return I
        
        YpxJa_col = YpxJa_above.apply(lambda x: Int(x,YpxJa_above.columns), axis = 1)
        YpxJp_col = YpxJp_above.apply(lambda x: Int(x,YpxJp_above.columns), axis = 1)
        
        tot_col = YpxJa_col + YpxJp_col
        I = np.trapz(tot_col.values[1:], tot_col.index[1:]) + tot_col.values[0]
        lat_int.append(I)
    
    
    lats_int = lat_int * np.cos(lats)
    Glob_prod = np.trapz(lats_int, lats) / 2
    phis_glob.append(Glob_prod)
    
  
