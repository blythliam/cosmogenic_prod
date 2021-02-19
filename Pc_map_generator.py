# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 10:31:24 2020

@author: lpb20
"""

#%%

####################################
##########Functions#################
####################################
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
    """
    This function creates the LIS for the given input energies depending on the 
    value of phi and whether desires for protons or alpha and heavier

    Parameters
    ----------
    E_vals : TYPE
        Energy Values for which the LIS needs to be created.
    phi : TYPE
        The local shielding potential in MeV.
    p_or_a : TYPE
        P=Proton / a = Alpha and heavy.

    Returns
    -------
    J : TYPE
        Returns J, the local insterellar spectrun for the energies, E provided.

    """
    Tr = 0.938
    if p_or_a == 'p':
        Z = 1
        A = 1
        phi = (Z/A) * phi
        T = E_vals
        Tphi = T + phi
        Gamma = np.array((Tr + Tphi) / Tr)
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

def Gen_PC(lats_points = 144, long_points = 192):
    """
    

    Parameters
    ----------
    lats_points : TYPE, optional
        The number of latitude points required. The default is 144.
    long_points : TYPE, optional
        The number of longitude points required. The default is 192.

    Returns
    -------
    Pc : TYPE
        A map of the shielding potential created by the earth's magnetic field, 
        given in MeV.

    """
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
    
    x0 = 399e3#((L1-g11*E) / 3*B0**2) * R0
    y0 = 351e3#((L2-h11*E) / 3*B0**2) * R0
    z0 = 221e3
    
    
    
    lats = np.linspace(0, pi, lats_points)
    longs = np.linspace(-pi, pi, long_points)
    
    long, lat = np.meshgrid(longs,lats)
    d = np.sqrt(x0**2 + y0**2 + z0**2)
    R2 = R0**2 + d - 2 * R0 * (x0 * np.sin(lat) * np.cos(long) + y0 * np.sin(lat) * np.sin(long) + z0*np.cos(lat)) 
    R = np.sqrt(np.abs(R2))
    cosphiG = (1/(R * B0)) * (g11*(x0 - R0 * np.sin(lat) * np.cos(long)) + h11 * (y0 - R0 * np.sin(lat) * np.sin(long)) + g10 * (z0 - R0 * np.cos(lat)))
    
    phiG = np.arccos(cosphiG)
    
    Pc = 1.9 * M * (R0 / R)**2 * np.sin(phiG)**4 / 1e9
    
    return Pc
    
    
    

    
    
#%%



import pandas as pd
import numpy as np
from math import pi
import time


plot = False
M= 7.8
Er = 0.938
glob_tot = False
spat_ocean = False

#%%
" Setting the global varibles for magnetic moment of earth and mass of proton"



M= 7.8
Er = 0.938

#%%







Cproddir = 'C:/Users/lpb20/OneDrive - Imperial College London/Documents/Odyssey/cosmogenic_prod/14CProduction'
prod_file_path= 'C:/Users/lpb20/OneDrive - Imperial College London/Documents/Odyssey/cosmogenic_prod/Output_Files'


if spat_ocean == True:
    Cproddir = '//home/lpb20/cosmogenic_prod/14CProduction'



#S_p = pd.read_excel('COproduction.xlsx', sheet_name='14C_p', skiprows=4, index_col=0)
S_p_load = np.loadtxt(Cproddir + '/COproduction_p.csv', delimiter= ',', skiprows=4, dtype = np.str)[1:,1:]
S_a_load = np.loadtxt(Cproddir + '/COproduction_a.csv', delimiter= ',', skiprows=4, dtype = np.str)[1:,1:]
S_p = np.zeros([S_p_load.shape[0],S_p_load.shape[1]])
S_a = np.zeros([S_a_load.shape[0],S_a_load.shape[1]])
for row in range(S_p.shape[0]):
    for col in range(S_p.shape[1]):
        if S_p_load[row, col] == 'N/A':
            S_p[row,col] = 0
        else:
            
            S_p[row,col] = float(S_p_load[row,col])
        if S_a_load[row, col] == 'N/A':
            S_a[row,col] = 0
        else:
            S_a[row,col] = float(S_a_load[row,col])
#%%

height_vals = np.loadtxt(Cproddir + '/COproduction_p.csv', delimiter=',', skiprows=5,usecols=[0])
energy_vals = np.loadtxt(Cproddir + '/energy_values.csv', delimiter=',')


phi_data = np.loadtxt(Cproddir + '/phi_over_time.txt',delimiter = '  ',skiprows = 9)

S_p = pd.DataFrame(S_p, columns=energy_vals,index = height_vals)
S_a = pd.DataFrame(S_a, columns=energy_vals,index = height_vals)



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
#    S.fillna(0,inplace = True)▼
Y_p = pi * Sp_smooth # Y is the smoothed version of S * pi 
Y_a = pi * Sa_smooth


#%%

import glob
no_of_months_done = len(glob.glob('C:/Users/lpb20/OneDrive - Imperial College London/Documents/Odyssey/cosmogenic_prod/Output_Files/*.txt'))
phis = (phi_data[1] / 1000)[no_of_months_done:]
Pcs = Gen_PC()
start_year = time.time()
CO_prod = np.zeros([10, 110, 144, 192])
for year, phi in enumerate(phis):
    print('----------------------------------------------------')
    print('The code is now going to iterate through year: ', + (year))
    print('----------------------------------------------------')
    J_p = LIS(Y_p.columns, phi, 'p') 
    J_a = LIS(Y_a.columns, phi, 'a')
    
    YpxJp = Y_p.apply(lambda x: x*J_p,axis = 1)
    YpxJa = Y_a.apply(lambda x: x*J_a,axis = 1)
    
    
    column_tot = np.zeros([144,192])
    month_slice = np.zeros([110,144,192])
    for i in range(144):
        
        if i%2 == 0 and i>0:
            end = time.time()
            print('There are '+ str(144 - i)+ ' left in the year '+ str(year))
            print('That interation took ' + str(end - start) + 'seconds')
        start = time.time()
        for j in range(192):


            Pc = Pcs[i,j]
            Ecp = 0.938 * (np.sqrt(1+((Pc)/(0.938))**2)-1)
            Eca = 0.938 * (np.sqrt(1 + ((2 * Pc)/(4 *0.938))**2) -1)
            
            YpxJa_above = YpxJa[YpxJa.columns[YpxJa.columns>Eca]]
            
            YpxJp_above = YpxJp[YpxJp.columns[YpxJp.columns>Ecp]]
            
            
            def Int(y,x):
                I = np.sum(np.nan_to_num(np.diff(y*x) / ((np.diff(np.log(y))/np.diff(np.log(x))) + 1),nan=0.0))
                return I
            
            YpxJa_col = YpxJa_above.apply(lambda x: Int(x[x>0],YpxJa_above.columns[x>0]), axis = 1)
            YpxJp_col = YpxJp_above.apply(lambda x: Int(x[x>0],YpxJp_above.columns[x>0]), axis = 1)
            
            tot_col = pd.Series(YpxJa_col.values + YpxJp_col.values, index = YpxJp_col.index)
            month_slice[:,i,j] = tot_col
            if glob_tot == True:
                I = np.trapz(tot_col.values[1:], tot_col.index[1:]) + tot_col.values[0]
                column_tot[i,j] = I
    for ii in range(110):
        X = month_slice[ii,:,:]
        np.savetxt(prod_file_path + '_' + str(year) + '_' + str(ii).zfill(3),X)          

    CO_prod[year, :,  : , : ] = month_slice
    end_year = time.time()
    
print('---------------------------------------------------')
print('\tOne Iteration took: ' + str((end_year-start_year)/60) + ' mins')
print('---------------------------------------------------')

      
  
