# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 11:03:55 2024

@author: elg
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import cartopy.feature as cf
import cartopy.crs as ccrs
import pygrib
import matplotlib.colors as mcolors


######### EVENT AM  Height and Wind Barbs
grbGFS = pygrib.open('gfs.0p25.2021121012.f000.grib2')

#Wind surface map with wind barbs
grbGFS.select(name='Geopotential Height')
grbGFS.select(name='U component of wind')
grbGFS.select(name='V component of wind')
grbGFS.select(name='Temperature')

# Read in data
gpm300 = grbGFS[115]; gpm = gpm300['values']
uwnd300 = grbGFS[121]; uwnd = uwnd300['values']
vwnd300 = grbGFS[122]; vwnd = vwnd300['values']
tmp300 = grbGFS[116]; tmp = tmp300['values']
lats, lons = gpm300.latlons()
lats, lons = uwnd300.latlons() 
lats, lons = vwnd300.latlons()
lats, lons = tmp300.latlons()

#Set up cartopy map
fig = plt.figure (figsize=(8,8))
proj=ccrs.LambertConformal(central_longitude=-96.,central_latitude=40.,standard_parallels=(40.,40.))
ax=plt.axes(projection=proj)
ax.set_extent([-125.,-70.,20.,60.])

#Add features
ax.add_feature(cf.LAND,color='white')
ax.add_feature(cf.OCEAN, color = 'white')
ax.add_feature(cf.COASTLINE,edgecolor='gray')
ax.add_feature(cf.STATES,edgecolor='gray')
ax.add_feature(cf.BORDERS,edgecolor='gray',linestyle='-')
ax.add_feature(cf.LAKES, alpha=0.5)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')

gl.top_labels = False
gl.left_labels = False

#Plot data
#Geopotential Height 
# Plot the 300mb geopotential heights as thick solid black contours
gpm_contour = plt.contour(lons, lats, gpm, 
                          levels=np.arange(np.min(gpm), np.max(gpm), 60),  # Adjust the levels as needed
                          colors='black', 
                          linewidths=2,  # Specify the line width
                          transform=ccrs.PlateCarree())


# Plot the temperature contour as shaded contours
tmp_contour = plt.contourf(lons, lats, tmp - 273, 
                           np.arange(np.min(tmp - 273), np.max(tmp - 273), 5), 
                           cmap=plt.cm.rainbow,  
                           transform=ccrs.PlateCarree())

cbar = plt.colorbar(tmp_contour)
cbar.set_label('Temperature (°C)')

plt.barbs(lons[::20,::20],lats[::20,::20],uwnd[::20,::20],vwnd[::20,::20],transform=ccrs.PlateCarree())

plt.title ('Temperature and Wind Barbs 300 mb Heights Dec 10, 2021 12 Z')
plt.savefig('FinalProjectmap1.png', dpi=300)
plt.show()
plt.close()

######### EVENT TIME  Height and Wind Barbs
grbGFS = pygrib.open('gfs.0p25.2021121100.f000.grib2')

#Wind surface map with wind barbs
grbGFS.select(name='Geopotential Height')
grbGFS.select(name='U component of wind')
grbGFS.select(name='V component of wind')
grbGFS.select(name='Temperature')

# Read in data
gpm300 = grbGFS[115]; gpm = gpm300['values']
uwnd300 = grbGFS[121]; uwnd = uwnd300['values']
vwnd300 = grbGFS[122]; vwnd = vwnd300['values']
tmp300 = grbGFS[116]; tmp = tmp300['values']
lats, lons = gpm300.latlons()
lats, lons = uwnd300.latlons() 
lats, lons = vwnd300.latlons()
lats, lons = tmp300.latlons()

#Set up cartopy map
fig = plt.figure (figsize=(8,8))
proj=ccrs.LambertConformal(central_longitude=-96.,central_latitude=40.,standard_parallels=(40.,40.))
ax=plt.axes(projection=proj)
ax.set_extent([-125.,-70.,20.,60.])

#Add features
ax.add_feature(cf.LAND,color='white')
ax.add_feature(cf.OCEAN, color = 'white')
ax.add_feature(cf.COASTLINE,edgecolor='gray')
ax.add_feature(cf.STATES,edgecolor='gray')
ax.add_feature(cf.BORDERS,edgecolor='gray',linestyle='-')
ax.add_feature(cf.LAKES, alpha=0.5)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')

gl.top_labels = False
gl.left_labels = False

#Plot data
#Geopotential Height 
# Plot the 300mb geopotential heights as thick solid black contours
gpm_contour = plt.contour(lons, lats, gpm, 
                          levels=np.arange(np.min(gpm), np.max(gpm), 60),  # Adjust the levels as needed
                          colors='black', 
                          linewidths=2,  # Specify the line width
                          transform=ccrs.PlateCarree())


# Plot the temperature contour as shaded contours
tmp_contour = plt.contourf(lons, lats, tmp - 273, 
                           np.arange(np.min(tmp - 273), np.max(tmp - 273), 5), 
                           cmap=plt.cm.rainbow,  
                           transform=ccrs.PlateCarree())

cbar = plt.colorbar(tmp_contour)
cbar.set_label('Temperature (°C)')

plt.barbs(lons[::20,::20],lats[::20,::20],uwnd[::20,::20],vwnd[::20,::20],transform=ccrs.PlateCarree())

plt.title ('Temperature and Wind Barbs 300 mb Heights Dec 11, 2021 03 Z')
plt.savefig('FinalProjectmap2.png', dpi=300)
plt.show()
plt.close()




