# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 14:56:45 2024

@author: Davis
"""

import cartopy # imports all of cartopy
import cartopy.feature as cf #imports the features class
import cartopy.crs as ccrs #imports the coordinate referencesystem class
import matplotlib.pyplot as plt #matplotlib is needed to plot the data
import pygrib # pygrib is needed to extract the data
import numpy as np # numpy is needed to do calculations on the data

December10grb = pygrib.open ('gfs.0p25.2021121012.f000.grib2') # Open the December 10th, 12z File
December11grb = pygrib.open ('gfs.0p25.2021121100.f000.grib2') # Open the December 11th, 0z File

# December 10th 300mb Winds

# Make a map of the US
fig = plt.figure
proj=ccrs.LambertConformal(central_longitude=-96.,central_latitude=30.,standard_parallels=(30.,30.))
ax=plt.axes(projection=proj)
ax.set_extent([-125.,-70.,20.,50.])
ax.add_feature(cf.LAND)
ax.add_feature(cf.OCEAN)
ax.add_feature(cf.COASTLINE)
ax.add_feature(cf.STATES)
ax.add_feature(cf.BORDERS)
ax.add_feature(cf.LAKES)
ax.add_feature(cf.STATES,edgecolor='grey')
ax.gridlines(draw_labels=True, linewidth=1, color='white', alpha=0.5, linestyle='--')

# Read the data
December10heights = December10grb.select(name='Geopotential height', level=300)[0]
December10u_wind = December10grb.select(name='U component of wind', level=300)[0]
December10v_wind = December10grb.select(name='V component of wind', level=300)[0]
lats, lons = December10heights.latlons()
December10wind_magnitude = np.sqrt(December10u_wind.values ** 2 + December10v_wind.values ** 2)

# Plot the data
cs_heights = ax.contour(lons, lats, December10heights.values / 10, transform=ccrs.PlateCarree(), levels=20, colors='black')
cs_wind = ax.contourf(lons, lats, December10wind_magnitude, transform=ccrs.PlateCarree(), levels=20, cmap='viridis')

# Add legends and title
cbar_wind = plt.colorbar(cs_wind, ax=ax, orientation='horizontal', shrink=0.5)
cbar_wind.set_label('Wind Magnitude (knots)')
plt.title('December 10th, 12z 300 mb Heights (dm) / Isotachs (knots)')

plt.show()


# December 10th 300mb Absolute Vorticity

# Read the Data
December10vorticity = December10grb.select(name="Absolute vorticity", level=300)[0]

# Make a map of the US
fig = plt.figure
proj=ccrs.LambertConformal(central_longitude=-96.,central_latitude=30.,standard_parallels=(30.,30.))
ax=plt.axes(projection=proj)
ax.set_extent([-125.,-70.,20.,50.])
ax.add_feature(cf.LAND)
ax.add_feature(cf.OCEAN)
ax.add_feature(cf.COASTLINE)
ax.add_feature(cf.STATES)
ax.add_feature(cf.BORDERS)
ax.add_feature(cf.LAKES)
ax.add_feature(cf.STATES,edgecolor='grey')
ax.gridlines(draw_labels=True, linewidth=1, color='white', alpha=0.5, linestyle='--')



# Plot the Data
cs_vorticity = ax.contourf(lons, lats, December10vorticity.values * 1e5, transform=ccrs.PlateCarree(), levels=2\0, cmap='YlOrBr')
cs_heights = ax.contour(lons, lats, December10heights.values / 10, transform=ccrs.PlateCarree(), levels=20, colors='black')

# Add legend and title
#plt.legend(*cs_vorticity.legend_elements(), title='10^-5 s^-1', loc='center left', bbox_to_anchor=(1, 0.5))
cbar_vort10 = plt.colorbar(cs_vorticity, ax=ax, orientation='horizontal', shrink=0.5)
plt.title('December 10th, 12z 500mb Heights (dm) / Absolute Vorticity (10^-5)')

plt.show()



# December 11th 300mb Winds

# Make a map of the US
fig = plt.figure
proj=ccrs.LambertConformal(central_longitude=-96.,central_latitude=30.,standard_parallels=(30.,30.))
ax=plt.axes(projection=proj)
ax.set_extent([-125.,-70.,20.,50.])
ax.add_feature(cf.LAND)
ax.add_feature(cf.OCEAN)
ax.add_feature(cf.COASTLINE)
ax.add_feature(cf.STATES)
ax.add_feature(cf.BORDERS)
ax.add_feature(cf.LAKES)
ax.add_feature(cf.STATES,edgecolor='grey')
ax.gridlines(draw_labels=True, linewidth=1, color='white', alpha=0.5, linestyle='--')

# Read the data
December11heights = December11grb.select(name='Geopotential height', level=300)[0]
December11u_wind = December11grb.select(name='U component of wind', level=300)[0]
December11v_wind = December11grb.select(name='V component of wind', level=300)[0]
lats, lons = December11heights.latlons()
December11wind_magnitude = np.sqrt(December11u_wind.values ** 2 + December11v_wind.values ** 2)

# Plot the data
cs_heights = ax.contour(lons, lats, December11heights.values / 10, transform=ccrs.PlateCarree(), levels=20, colors='black')
cs_wind = ax.contourf(lons, lats, December11wind_magnitude, transform=ccrs.PlateCarree(), levels=20, cmap='viridis')

# Add legends and title
cbar_wind = plt.colorbar(cs_wind, ax=ax, orientation='horizontal', shrink=0.5)
cbar_wind.set_label('Wind Magnitude (knots)')
plt.title('December 11th, 0z 300 mb Heights (dm) / Isotachs (knots)')

plt.show()


# December 11th 300mb Absolute Vorticity

# Read the Data
December11vorticity = December11grb.select(name="Absolute vorticity", level=300)[0]

# Make a map of the US
fig = plt.figure
proj=ccrs.LambertConformal(central_longitude=-96.,central_latitude=30.,standard_parallels=(30.,30.))
ax=plt.axes(projection=proj)
ax.set_extent([-125.,-70.,20.,50.])
ax.add_feature(cf.LAND)
ax.add_feature(cf.OCEAN)
ax.add_feature(cf.COASTLINE)
ax.add_feature(cf.STATES)
ax.add_feature(cf.BORDERS)
ax.add_feature(cf.LAKES)
ax.add_feature(cf.STATES,edgecolor='grey')
ax.gridlines(draw_labels=True, linewidth=1, color='white', alpha=0.5, linestyle='--')

# Plot the Data
cs_vorticity = ax.contourf(lons, lats, December11vorticity.values * 1e5, transform=ccrs.PlateCarree(), levels=10, cmap='coolwarm')
cs_heights = ax.contour(lons, lats, December11heights.values / 10, transform=ccrs.PlateCarree(), levels=10, colors='black')

# Add legend and title
cbar_vort11 = plt.colorbar(cs_vorticity, ax=ax, orientation='horizontal', shrink=0.5)
plt.title('December 11th, 12z 500mb Heights (dm) / Absolute Vorticity (10^-5)')

plt.show()