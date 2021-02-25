# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 14:11:47 2021

@author: lpb20
"""

" Converts the txt files representing the production of CO to a netCDF file"

import netCDF4 as nc
import numpy as np
#%%
" 14C production in atmosphere "







#%%
" Loading data and converting into netCDF file format"




df = nc.Dataset('CO_prod.nc', 'w',format = 'NETCDF4' )

times = df.createDimension('time', None)
lats_nc = df.createDimension('latitude', 144)
longs_nc = df.createDimension('longitude', 192)
heights_nc = longs_nc = df.createDimension('height', 110)
times = df.createVariable('time','f4',('time',))
lats_nc = df.createVariable('latitude', 'f4', ('latitude',))
longs_nc = df.createVariable('longitude', 'f4', ('longitude',))
heights_nc = df.createVariable('height', 'f4', ('height',))
value = df.createVariable('value','f4',('time','height','latitude','longitude'))
value.units = 'atoms cm-2'
lats_nc[:]=list(np.linspace(-90,90,144))
longs_nc[:] = list(np.linspace(-180,180,192))
times[:]  = [1,2]
heights_nc[:] = np.linspace(0,1000,110)
for year in range(no_of_years):
    for height in range(110): 
        value[year, height, : , : ] = CO_prod[year, height, : , : ]


df.close()


#%%
" Converting metadata to values required for input into UKCA"























