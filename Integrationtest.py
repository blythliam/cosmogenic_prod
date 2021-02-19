# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 08:39:54 2020

@author: blyth
"""
#Importing modules
import pandas as pd
import cairo
import matplotlib.pyplot as plt
#%matplotlib qt
import numpy as np
from math import pi

#Deciding whether to print plots or not
plot = False
M= 7.8
#%%
#Reading in data
p_a = 'a'
data = pd.read_excel('COproduction.xlsx', sheet_name='10Be_'+p_a, skiprows=4, index_col=0)
#Plotting Source function as function of altidude for specific energy
#data = pd.read_csv('Clproduction.csv', skiprows=4, index_col=0).iloc[:-1]
data.columns = data.columns.astype(float)
if plot:
    ploty = data[0.1]

    plt.loglog(list(data.index[1:]),ploty[1:])
    axes_1 = plt.axes()
    axes_1.set_ylim([1e-9,1e-4])
    axes_1.set_xlim([1,250])
    plt.show()

#Y = pi * data 

#%%
#Using trapz to integrate source function along height to find columnar production
integrate_test = data.apply(lambda x: np.trapz(x[1:][x[1:].notnull()],
                                               data.index[1:len(x[1:][x[1:].notnull()])+1])+max(x))
                                               
                                              
#%%
#Plotting columnar production for different energy

if plot:
    plt.figure(2,figsize=(5,15))    
    axes = plt.axes()
    axes.tick_params(labelright = True, right = True)
    plt.loglog(list(integrate_test.index),pi * np.array(integrate_test.values))
#    axess = plt.axes()
#    axess.set_ylim([0.001,100])
#    axess.set_xlim([0.1,100])
    plt.ylim([1e-3,100])
    plt.xlim([0.02,100])
    plt.xlabel('E(GeV)')
    plt.ylabel('Columnar Production')
    plt.title('Columnar production for 14C (protons)')
    
#%%

Yi = data.iloc[-1]
if plot:
    plt.loglog(Yi.index, Yi.values)
    plt.xlabel('E(Gev')
    plt.ylabel('$Y_i(E,1000)$')
    plt.title('The yield function at an atmsopheric depth of'+
                  ' 1000g/cm2 (surface)')
    
Y = data * pi
#YixJ = Y * J
#%%
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
#T = np.linspace(0,100,10000000)

#%% Setting up the modulated incoming spectra
    T = np.array(data.columns) #Using the energies from data not own energies
#    T = np.linspace(0,1000,10000)
    Tphi = T + phi_o
    
    PTphi = (Tphi*(Tphi+2*T_r))**0.5
    if p_a == 'p':
        JLIS = (1.9e4 * PTphi**-2.78) / (1 + 0.4866 * PTphi**-2.51)
    if p_a == 'a':
        JLIS = 0.05 * (1.9e4 * PTphi**-2.78) / (1 + 0.4866 * PTphi**-2.51)
    
    J = JLIS * T * (T + 2 * T_r) / (T + phi_o) / (T + phi_o + 2*T_r)/10000
    
    
    
    if plot:
        plt.figure(3)
    #    plt.loglog(T,JM)
        plt.loglog(T,J)
        plt.xlabel('E(GeV)')
        plt.ylabel('Intensity')
        plt.title('Plot of LIS at a specific time [Burger et al (2000)]')
        plt.ylim([0.01,3500])
        plt.xlim([0.1,100])
        
    #Finding the function YixJ which can be integrated to find Q.
    Y = (pi * data)
    Y.fillna(0,inplace = True)
    # Y2 = np.zeros([110,10000])
    # array2 = np.linspace(0,1000,10000)
    # for i in range(110):
    #     m, b = np.polyfit(Y.iloc[i].index,Y.iloc[i].values,1)
    #     Y2[i,:] = m * array2 + b
    
    YixJ = Y * J

    YixJ = pd.DataFrame(YixJ)
    
    


    

    #%%
    #Creating the grid of lat-long for plotting and the cutt off energies for 
    #the lower bound of the integral
    
#    latlon = np.linspace(-pi/2,pi/2,180)
    
#    Pc = 1.9 * M * (np.cos(latlon))**4
    
#    Ec = 0.938 * (np.sqrt(1+(Pc/0.938)**2)-1)
#    if plot:
#        plt.figure()
#        plt.plot(latlon*180/pi, Ec)
#        plt.xlabel('Latitide')
#        plt.ylabel('$E_{c}(GeV)$')
    
    
    
    #%%
    
    
    # finalint = {}
    # for lat in latlon:
        
    #     Pc = 1.9 * M * (np.cos(lat))**4
    #     Ec = T_r * (np.sqrt(1+((Pc/0.938)**2))-1)
    #     YixJ_lat = YixJ[YixJ.index>Ec]
    #     YixJ_lat = YixJ_lat[YixJ_lat.notnull()]
    #     finalint[lat] = np.trapz(YixJ_lat.values,YixJ_lat.index)
    #     finalint[-lat] = np.trapz(YixJ_lat.values,YixJ_lat.index)
            
            
    # final_int = pd.Series(finalint).sort_index()
    # final_int.index = final_int.index *180/pi
    # if plot:
    #     plt.plot(final_int.index,final_int.values)
    
    
    
    
    #%%
        
        
    # from mpl_toolkits.basemap import Basemap
    # lons, lats = np.meshgrid(lon,lat)  
    # if plot:
    #     mp = Basemap()
    #     x, y = mp(lon, lat)
    #     mp.drawcoastlines()
    #     c_scheme = mp.pcolor(x, y, produc, cmap = 'jet') 
    #     cbar = mp.colorbar(c_scheme)
    
        
    #%%
    
    
    
    
    
    #%%
    # Producing a YixJ for all lats and longs
    
    num_points = 180   
    lats = np.transpose(np.linspace(-pi/2,pi/2,num_points))
    
    
    Q_plane = np.zeros((num_points,110,23))
    for i in range(num_points):
        Q_plane[i,:,:] = YixJ
        
    lat_int = np.zeros([num_points])
    
    for i,lat in enumerate(lats):
        if p_a == 'a':
            Emin = 0.08
        if p_a == 'p':
            Emin = 0.02
            
        Pc = 1.9 * M * (np.cos(lat))**4
        Ec = 0.938 * (np.sqrt(1+(Pc/0.938)**2)-1)
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
        lat_int[i] =  np.trapz(height_int[1:].values,height_int[1:].index)+height_int[:1].values    
            
        

# lat_plane.insert(0,Ec,np.zeros([110,1]))
    #       height_int = lat_plane.apply(lambda x: np.trapz(x[x.notnull()],
    #                               x[x.notnull()].index),axis = 1)
        
        
    
    final_int_val = (np.trapz(lat_int * np.cos(lats),lats) * 2 * pi) / (4 * pi) 
    prod.append(final_int_val)

#%%  Unsing approach from Asvertari et al 2017


# columnar_prod_low = Y[1:].apply(lambda x: np.trapz(x[x.notnull()],x[x.notnull()].index))

# Yc = columnar_prod_low + Y[:1]

# #plt.loglog(np.array(columnar_prod.columns).reshape((23,1)),
# #           np.array(columnar_prod.values).reshape((23,1))/pi)

# YcxJ = Yc * J 

# lambda_ = np.linspace(0 , pi/2 , 180)

# Q_plane_2 = np.zeros((180,23))
# for j in range(180):
#     Q_plane_2[j,:] = YcxJ
# lat_val =[]
# for i,lat in enumerate(lambda_):
#     lat_slice = pd.Series(Q_plane_2[i,:]) * np.cos(lat)
#     lat_slice.index = YcxJ.columns
#     Pc = 1.9 * M * (np.cos(lat))**4
#     Ec = 0.938 * (np.sqrt(1+(Pc/0.938)**2)-1)
    
#     lat_val.append(np.trapz(lat_slice[lat_slice.index>Ec].values,lat_slice[lat_slice.index>Ec].index))
    
# plt.plot(phi_os,M_12,'-.', phi_os, M_78, phi_os, M_6,'--')
# plt.legend(['M = 12', 'M = 7.8', 'M = 6'])
# plt.ylim([0,3.4])













    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


