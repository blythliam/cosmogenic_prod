# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 10:24:05 2020

@author: lpb20
"""

import pandas as pd
import cairo
import matplotlib.pyplot as plt
#%matplotlib qt
import numpy as np
from math import pi

#Deciding whether to print plots or not
input_mol = '14C'
plot = True
for M in [6,7.8,12]:
    
#Reading in data
    for p_a in ['p','a']:
        data = pd.read_excel('COproduction.xlsx', sheet_name=input_mol+'_'+p_a, skiprows=4, index_col=0)
        #Plotting Source function as function of altidude for specific energy
        #data = pd.read_csv('Clproduction.csv', skiprows=4, index_col=0).iloc[:-1]
        data.columns = data.columns.astype(float)
        
        
        #Creating Local Interstellar spectrum 
        if p_a == 'p':
            Z = 1
            A = 1
        if p_a == 'a':
            Z = 2
            A = 4
        T_r = 0.938
        
        phi_os = np.linspace(0,2,10)
        prod = []
        for phi_o in phi_os:
            phi_o = phi_o * (Z/A)
            print(phi_o)
            T = np.array(data.columns) #Using the energies from data not own energies
        #    T = np.linspace(0,1000,10000)
            Tphi = T + phi_o
            
            PTphi = (Tphi*(Tphi+2*T_r))**0.5
            if p_a == 'p':
                JLIS = (1.9e4 * PTphi**-2.78) / (1 + 0.4866 * PTphi**-2.51)
            if p_a == 'a':
                JLIS = 0.05 * (1.9e4 * PTphi**-2.78) / (1 + 0.4866 * PTphi**-2.51)
            
            J = JLIS * T * (T + 2 * T_r) / (T + phi_o) / (T + phi_o + 2*T_r)/10000
            
            
            #Finding the function YixJ which can be integrated to find Q.
            Y = (pi * data)
            Y.fillna(0,inplace = True)
            YixJ = Y * J
            YixJ = pd.DataFrame(YixJ)
            
        
            # Producing a YixJ for all lats and longs
            
            num_points = 180
            lats = np.transpose(np.linspace(-pi/2 , pi/2 , num_points))
            
            
            Q_plane = np.zeros((num_points,110,23))
            for i in range(num_points):
                Q_plane[i,:,:] = YixJ
                
            lat_int = np.zeros([num_points])
            height_ints = np.zeros([110,num_points])
            for i,lat in enumerate(lats):
                if p_a == 'a':
                    Emin = 0.08
                if p_a == 'p':
                    Emin = 0.02
                    
                Pc = 1.9 * M * (np.cos(lat))**4
                Ec = 0.938 * (np.sqrt(1+((Pc)/(0.938))**2)-1)
                lat_plane = pd.DataFrame(Q_plane[i,:,:])
                if p_a == 'a':
                    lat_plane.columns =  YixJ.columns * 4
                elif p_a == 'p':
                    lat_plane.columns =  YixJ.columns
                lat_plane.index = YixJ.index
                if Ec > Emin:
                    bot_col = lat_plane[lat_plane.columns[lat_plane.columns>Ec]].iloc[:,-0]
                    top_col = lat_plane[lat_plane.columns[lat_plane.columns<Ec]].iloc[:,-1]
                lat_plane = lat_plane[lat_plane.columns[lat_plane.columns > Ec]] 
                if Ec>Emin:
                    insert_col = abs(bot_col - top_col)/2
                    lat_plane.insert(0,Ec,insert_col)
                     
                height_int = lat_plane.apply(lambda x: np.trapz(x,x.index),axis = 1)
                height_ints[:,i] = height_int.values
                lat_int[i] =  np.trapz(height_int[1:].values,height_int[1:].index)+height_int[:1].values    
                    
            
            final_int_val = (np.trapz(lat_int * np.cos(lats),lats) * 2 * pi) / (4 * pi) 
            
            prod.append(final_int_val)
            
        if p_a == 'p':
            prot_prod = prod
            prot_height = height_ints
            lat_int_prot = lat_int
        if p_a == 'a':
            prot_alpha = prod
            alph_height = height_ints
            lat_int_alpha = lat_int
            
    if plot:
        prod_tot = np.array(prot_prod) + np.array(prot_alpha)
        if M == 6:
            plt.plot(phi_os, prod_tot,'--')
        elif M == 7.8:
            plt.plot(phi_os, prod_tot)
        elif M == 12:
            plt.plot(phi_os, prod_tot,'.-')





plt.ylim([0.0, 3.4])
plt.xlim([0,2.0])
plt.legend(['M = 6, M = 7.8, M = 12'])
plt.ylabel('Global Production (at/cm$^2$/sec)')
plt.xlabel('Modulation Potential $\phi$ (GeV)')
plt.title('Global Production of '+input_mol)

    
#    np.append(prod, final_int_val) 

#%%






if plot:
    height = np.array(data.index)
    height[0] = 5
    height_func = prot_height + alph_height
    height_func = height_func / prod_tot
    height_func = height_func / 1e-3
    levels = [0.2, 0.4, 0.6, 0.8, 2, 4, 6, 8, 10]
    cp = plt.contour(lats * (180/pi), height / 1000, height_func, levels, colors = 'k')
    plt.clabel(cp, inline =True)
    plt.clim(0,10)
    plt.gca().invert_yaxis()
    plt.yscale('log')
    plt.xlabel('Latitude')
    plt.ylabel('Atmospheric depth / 1000')
    

#plt.colorbar()
#if plot:
#    plt.contour
#%%

lat_plane = np.zeros([180,5]) 
longs = np.linspace(-180,180,5)
lats = np.linspace(-90,90,180)
for i in range(5):
    lat_plane[:,i] = ((lat_int_prot +lat_int_alpha))/2
    



levels = [0.4, 0.8, 1.2, 1.6, 2]
cp = plt.contour(lat_plane, levels, colors = 'k')
plt.clabel(cp, inline =True)
plt.clim(0,10)
plt.gca().invert_yaxis()



























   