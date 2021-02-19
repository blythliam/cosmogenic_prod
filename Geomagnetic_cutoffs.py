# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 11:08:00 2021

@author: lpb20
"""



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import pi

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

L0 = 2* g10 * g20 + (np.sqrt(3) * (g11*g21 + h11*h21))
L1 = (-1 * (g11 * g20)) + (np.sqrt(3) * (g10*g21 + g11*g22 + h11*h22))
L2 = -1 * h11 * g20 + (np.sqrt(3) * (g10*h21 - h11*g22 + g11*h22))

E = ((L0*g10 + L1 * g11  + L2 * h11) / 4 * B0**2)

x0 = -399e3#((L1-g11*E) / 3*B0**2) * R0
y0 = 351e3#((L2-h11*E) / 3*B0**2) * R0
z0 = 221e3#((L0-g10*E) / 3*B0**2) * R0



lats = np.linspace(0, pi, 144)
longs = np.linspace(-pi, pi, 192)

long, lat = np.meshgrid(longs,lats)
d = np.sqrt(x0**2 + y0**2 + z0**2)
R2 = R0**2 + d - 2 * R0 * (x0 * np.sin(lat) * np.cos(long) + y0 * np.sin(lat) * np.sin(long) + z0*np.cos(lat)) 
R = np.sqrt(np.abs(R2))
cosphiG = (1/(R * B0)) * (g11*(x0 - R0 * np.sin(lat) * np.cos(long)) + h11 * (y0 - R0 * np.sin(lat) * np.sin(long)) + g10 * (z0 - R0 * np.cos(lat)))

phiG = np.arccos(cosphiG)
Geophi = pi/2 - phiG
Pc = (1.9 * M * (R0 / R)**2 * np.cos(Geophi)**4) / 1e9
Pc_approx = 1.9 * M * (np.cos(lat - pi/2)) ** 4 / 1e9
Pc = np.nan_to_num(Pc)
plt.figure()
plt.imshow(Pc, extent = [-180, 180, -90, 90])
cnorm = ['']
cbar = plt.colorbar()
cbar.set_label('P$_c$ (Gev)')
plt.title('Values of P$_c$ for the Eccentric Dipole approximate')

#%%

fig = plt.figure()
ax = fig.gca(projection = '3d')

surf = ax.plot_surface(lat, long, Pc)

#%%

Pc = np.roll(Pc, 45, axis = 1)



values = np.zeros([144*192,1])
lats_df = np.zeros([144*192,1])
longs_df = np.zeros([144*192,1])

for i in range(144):
    for j in range(192):
        values[(192 * i) +j] = Pc[i,j]
        lats_df[(192*i + j)] = (-1 * (lats[i] * 180/pi)) + 90
        longs_df[192*i + j] = longs[j] * 180/pi

export_file = pd.DataFrame(np.concatenate((lats_df,longs_df,values), axis = 1))
export_file.columns = ['lat','lon','Pc']

export_file.to_csv('Pc.csv', index = False)


#%%


plt.figure()
plt.imshow(Pc_approx, extent = [-180, 180, -90, 90])
cbar = plt.colorbar()
cbar.set_label('P$_c$ (Gev)')
plt.title('Values of P$_c$ for the Geocentric Dipole')



#%%

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
# lon_0 is central longitude of projection.
# resolution = 'c' means use crude resolution coastlines.
m = Basemap(projection='mbtfpq',lon_0=0,resolution='c')
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='white')
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,360.,60.))
m.drawmapboundary(fill_color='white')
plt.title("McBryde-Thomas Flat Polar Quartic Projection")


#plt.colorbar()





