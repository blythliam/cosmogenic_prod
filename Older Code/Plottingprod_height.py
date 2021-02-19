# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 16:38:33 2020

@author: lpb20
"""

# Trying to recretae similar plots to those in Pavlov et al 2017

import pandas as pd
import numpy as np
import cairo

import matplotlib.pyplot as plt
#%matplotlib qt
plot = True
from math import pi


data = pd.read_excel('COproduction.xlsx', sheet_name='14C_p', skiprows=4, index_col=0)

integrate_test = data.apply(lambda x: np.trapz(x[1:][x[1:].notnull()],
                                               data.index[1:len(x[1:][x[1:].notnull()])+1])+max(x))
   

if plot:
    plt.figure(2,figsize=(5,15))    
    axes = plt.axes()
    axes.tick_params(labelright = True, right = True)
    plt.loglog(list(integrate_test.index),pi * np.array(integrate_test.values))
    axes = plt.axes()
    axes.set_ylim([1e-3,1e2])
    axes.set_xlim([0.01,100])
    plt.xlabel('E(GeV)')
    plt.ylabel('Columnar Production')
    plt.title('Columnar production for 14C (proton)')
    plt.show()
    
    
R0 = 14.9
lamb = 79.36
R = R0 * (np.cos(lamb))**4

True_data = np.transpose(data.columns > R)


#data_up = data.where(data.columns > R)
















    

    
    
