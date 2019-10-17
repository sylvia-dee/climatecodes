	
# Calculate Climate Velocities For Gulf of Mexico.

import numpy as np
import numpy.ma as ma

import string
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import sys, os, cdtime, cdutil, cdms2, MV2, time, datetime
from math import exp
import scipy.io as sio


cd /rdf/sd75/sylvia/GULF_CESM/
T1=np.load('SST_ENSMEAN_2006-2080.npy')
T2=np.load('SST_ENSMEAN_2080-2100.npy')

LAT=np.load('POP2_lats.npy')
LON=np.load('POP2_lon.npy')

LON=LON[0,:]
LAT=LAT[:,0]
#======================================================================

#======================================================================

#======================================================================
# SAVE AS MEAN FROM 2006-2026, 2080-2100
#======================================================================

TEMP_1=np.mean(T1[0:240,:,:],axis=0)
TEMP_2=np.mean(T2,axis=0)

#======================================================================
# 1. TEMP RATE OVER TIME
#======================================================================
T1_mask=np.ma.masked_outside(TEMP_1,-10,400.)
T2_mask=np.ma.masked_outside(TEMP_2,-10,400.)

# CALCULATE RATE OF CLIMATE CHANGE, C/YR

T_DIFF1=np.mean(T1[228:240,:,:],axis=0)-np.mean(T1[0:12,:,:],axis=0)
T_DIFF2=np.mean(T2[-12::,:,:],axis=0)-np.mean(T2[0:12,:,:],axis=0)

# use centerpoint year to establish change over number of years [2006-2026], [2080-2100]
nyears=20.

Rate_time1=T_DIFF1/nyears
Rate_time2=T_DIFF2/nyears


Rate_time1=np.ma.masked_less(Rate_time1,0.001)
Rate_time2=np.ma.masked_less(Rate_time2,0.001)

#======================================================================
# 2. CALCULATE CHANGE IN TEMPERATURE S-N, W-E
#======================================================================

# Climate velocity is calculated by dividing the rate of climate change 
# by the rate of spatial climate variability to obtain a speed at which 
# species must migrate over the surface of the earth to maintain constant climate conditions

# 3.  Calculate the "centered" difference for T at each grid box, in both the east-west and north-south, direction, e.g.:

# T_diff_ew = R(x+1)-R(x-1)   where x is the east-west direction and "+1" means the grid box to the east and "-1" means the grid box to the west.

TC_diff_ew=np.zeros((T1_mask.shape))
TE_diff_ew=np.zeros((T1_mask.shape))

for i in range(383):
	for j in range(0,319):
		TC_diff_ew[i,j]=T1_mask[i,j+1]-T1_mask[i,j-1]
		TE_diff_ew[i,j]=T2_mask[i,j+1]-T2_mask[i,j-1]


# T_diff_ns  = R(y+1)-R(y-1)   where y is the north-south direction and "+1" means the grid box to the north and "-1" means the grid box to the south.

TC_diff_ns=np.zeros((T1_mask.shape))
TE_diff_ns=np.zeros((T1_mask.shape))


for i in range(383):
	for j in range(0,319):
		TC_diff_ns[i,j]=T1_mask[i+1,j]-T1_mask[i-1,j]
		TE_diff_ns[i,j]=T2_mask[i+1,j]-T2_mask[i-1,j]

# SAME AFTER THIS.
#======================================================================


# 4.  Calculate the change (difference) in the actual distance between grid points.  To do that, first just calculate the change in latitude and longitude:

# dlat = lat(y) - lat(y-1)  
# dlon = lon(x)-lon(x-1)

dlat=np.diff(LAT)
dlon=np.diff(LON)

# then convert into meters:

# dy = dlat*111.3 km * 1000 m/km
# dx = dlon*111.3 km * 1000 m/km * cos(lat)

# NOTE: np.cos takes amount in radians--have to convert lat to radians

rad_lat=LAT*np.pi/180.
rad_lat=np.float32(rad_lat)

dy=dlat[0]*111.3 # m
dy=np.float32(dy)

dx=np.zeros((383,319))
for i in range(319):
	for j in range(383):
		dx[j,i]=dlon[i]*111.3*np.cos(rad_lat[j])

# Basically the north-south distance (dy) is constant at every grid point, but the East-west distance (dx) 
# varies as a function of the cosine of the latitude.  
# In other words, dx should be an array that varies with latitude, while dy should just be a number.

# 5.  Calculate the spatial derivatives of the isotope ratio, i.e.:

TC_ew_deriv = TC_diff_ew[0:-1,0:-1]/(2.*dx) # UNITS: delta C / delta x (km)
TE_ew_deriv = TE_diff_ew[0:-1,0:-1]/(2.*dx)

TC_ns_deriv = TC_diff_ns[0:-1,0:-1]/(2.*dy)
TE_ns_deriv = TE_diff_ns[0:-1,0:-1]/(2.*dy)



	# x_gradient <- raster::focal(WE_gradient, w = f, nrow=3, ncol=3, pad = TRUE,
	# 											 fun = function(x, ...) {
	# 											 	#browser()
	# 											 	wt <- c(1, 2, 1, 2, 0, 2, 1, 2, 1)
	# 											 	weighted.mean(x, wt, na.rm = TRUE)
	# 											 })
	# y_gradient <- raster::focal(NS_gradient, w = f, nrow = 3, ncol = 3, pad = TRUE,
	# 											 fun = function(x, ...) {
	# 											 	wt <- c(1, 2, 1, 2, 0, 2, 1, 2, 1)
	# 											 	weighted.mean(x, wt, na.rm = TRUE)
	# 											 })


# Get magnitude of resultant vector

TC_MAG=np.sqrt(np.power(TC_ew_deriv,2) + np.power(TC_ns_deriv,2))
TE_MAG=np.sqrt(np.power(TE_ew_deriv,2) + np.power(TE_ns_deriv,2))

	# magnitude <- sqrt(x_gradient^2 + y_gradient^2)

	# # Get angle of resultant vector
	# angle <- raster::atan2(x_gradient, y_gradient) * 180 / pi


	# get_vocc <- function(hadsst_raster, years = 1969:2009, truncate = FALSE) {

	# linear_change <- get_sst_linear_change(hadsst_raster, years)

	# average_sst <- get_average_sst(hadsst_raster, years)

	# WE_gradient <- get_WE_diffs(average_sst)
	# NS_gradient <- get_NS_diffs(average_sst)

	# spatial_gradient <- get_spatial_gradient(NS_gradient, WE_gradient)

	# velocity <- linear_change / raster::subset(spatial_gradient, 'spatial_gradient')


#======================================================================
# CALCULTE CLIMATE VELOCITY:
# (C yr-1 / C km-1 = km / yr )

CLIM_VEL1=Rate_time1[0:-1,0:-1]/TC_MAG
CLIM_VEL2=Rate_time2[0:-1,0:-1]/TE_MAG
#======================================================================

# MAP,lons2 =m.shiftdata(SST_DIFF,Slon, lon_0=0)
from mpl_toolkits.basemap import Basemap, shiftgrid
slevels=np.arange(0,200,10)
tlevels=np.arange(0,0.05,0.005)
mlevels=np.arange(0,0.015,0.001)

fig,ax=plt.subplots()
#======================================================================

plt.subplot(3,2,1)
#m = Basemap(projection='robin',lon_0=0,resolution='c')
m = Basemap(projection='cyl', llcrnrlon=260, llcrnrlat=7, urcrnrlon=295, urcrnrlat=33)
# m = Basemap(width=6000000,height=4000000,projection='lcc',
#              resolution=None,lat_1=27.,lat_2=30.,lat_0=25.,lon_0=-90.)
# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(LON[0:-1], LAT[0:-1]))
#m.drawcoastlines()
# draw filled contours.
#ax1 = m.contourf(x,y,VAR_DIFF_mask,50,cmap=plt.cm.plasma_r,levels=slevels)
ax1 = m.contourf(x,y,CLIM_VEL1,50,cmap=plt.cm.magma_r,levels=slevels)

ax2 = m.contour(x,y,CLIM_VEL1,50,levels=slevels,colors='k')
plt.clabel(ax2, inline=1, fontsize=10)
cbar = m.colorbar(ax1,location='right',pad="2%")
#m.fillcontinents(color='coral',lake_color='aqua')
cbar.set_label(r'VELOCITY, [km/yr]')
fig.patch.set_visible(False)
ax.axis('off')
plt.title("A. CLIMATE VELOCITY, [2006-2026] [km/yr]")
#======================================================================
plt.subplot(3,2,2)
#m = Basemap(projection='robin',lon_0=0,resolution='c')
m = Basemap(projection='cyl', llcrnrlon=260, llcrnrlat=7, urcrnrlon=295, urcrnrlat=33)
# m = Basemap(width=6000000,height=4000000,projection='lcc',
#              resolution=None,lat_1=27.,lat_2=30.,lat_0=25.,lon_0=-90.)
# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(LON[0:-1], LAT[0:-1]))
#m.drawcoastlines()
# draw filled contours.
#ax1 = m.contourf(x,y,VAR_DIFF_mask,50,cmap=plt.cm.plasma_r,levels=slevels)
ax1 = m.contourf(x,y,CLIM_VEL2,50,cmap=plt.cm.magma_r,levels=slevels)

ax2 = m.contour(x,y,CLIM_VEL2,50,levels=slevels,colors='k')
plt.clabel(ax2, inline=1, fontsize=10)
cbar = m.colorbar(ax1,location='right',pad="2%")
#m.fillcontinents(color='coral',lake_color='aqua')
cbar.set_label(r'VELOCITY, [km/yr]')
fig.patch.set_visible(False)
ax.axis('off')
plt.title("B. CLIMATE VELOCITY, [2080-2100] [km/yr]")

#======================================================================
plt.subplot(3,2,3)
#m = Basemap(projection='robin',lon_0=0,resolution='c')
m = Basemap(projection='cyl', llcrnrlon=260, llcrnrlat=7, urcrnrlon=295, urcrnrlat=33)
# m = Basemap(width=6000000,height=4000000,projection='lcc',
#              resolution=None,lat_1=27.,lat_2=30.,lat_0=25.,lon_0=-90.)
# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(LON[0:-1], LAT[0:-1]))
#m.drawcoastlines()
# draw filled contours.
#ax1 = m.contourf(x,y,VAR_DIFF_mask,50,cmap=plt.cm.plasma_r,levels=slevels)
ax1 = m.contourf(x,y,Rate_time1[0:-1,0:-1],50,cmap=plt.cm.Reds,levels=tlevels)

ax2 = m.contour(x,y,Rate_time1[0:-1,0:-1],50,levels=tlevels,colors='k')
plt.clabel(ax2, inline=1, fontsize=10)
cbar = m.colorbar(ax1,location='right',pad="2%")
#m.fillcontinents(color='coral',lake_color='aqua')
cbar.set_label(r'Temperature / Time [$^\circ$C/yr]')
fig.patch.set_visible(False)
ax.axis('off')
plt.title("C. SST Change Rate, [2006-2026] [$^\circ$C/yr]")
#======================================================================
plt.subplot(3,2,4)
#m = Basemap(projection='robin',lon_0=0,resolution='c')
m = Basemap(projection='cyl', llcrnrlon=260, llcrnrlat=7, urcrnrlon=295, urcrnrlat=33)
# m = Basemap(width=6000000,height=4000000,projection='lcc',
#              resolution=None,lat_1=27.,lat_2=30.,lat_0=25.,lon_0=-90.)
# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(LON[0:-1], LAT[0:-1]))
#m.drawcoastlines()
# draw filled contours.
#ax1 = m.contourf(x,y,VAR_DIFF_mask,50,cmap=plt.cm.plasma_r,levels=slevels)
ax1 = m.contourf(x,y,Rate_time2[0:-1,0:-1],50,cmap=plt.cm.Reds,levels=tlevels)

ax2 = m.contour(x,y,Rate_time2[0:-1,0:-1],50,levels=tlevels,colors='k')
plt.clabel(ax2, inline=1, fontsize=10)
cbar = m.colorbar(ax1,location='right',pad="2%")
#m.fillcontinents(color='coral',lake_color='aqua')
cbar.set_label(r'Temperature / Time [$^\circ$C/yr]')
fig.patch.set_visible(False)
ax.axis('off')
plt.title("D. SST Change Rate, [2080-2100] [$^\circ$C/yr]")
#======================================================================

plt.subplot(3,2,5)
#m = Basemap(projection='robin',lon_0=0,resolution='c')
m = Basemap(projection='cyl', llcrnrlon=260, llcrnrlat=7, urcrnrlon=295, urcrnrlat=33)
# m = Basemap(width=6000000,height=4000000,projection='lcc',
#              resolution=None,lat_1=27.,lat_2=30.,lat_0=25.,lon_0=-90.)
# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(LON[0:-1], LAT[0:-1]))
#m.drawcoastlines()
# draw filled contours.
#ax1 = m.contourf(x,y,VAR_DIFF_mask,50,cmap=plt.cm.plasma_r,levels=slevels)
ax1 = m.contourf(x,y,TC_MAG,50,cmap=plt.cm.viridis_r,levels=mlevels)

ax2 = m.contour(x,y,TC_MAG,50,levels=mlevels,colors='k')
plt.clabel(ax2, inline=1, fontsize=10)
cbar = m.colorbar(ax1,location='right',pad="2%")
#m.fillcontinents(color='coral',lake_color='aqua')
cbar.set_label(r'Temperature / km, [$^\circ$C/km]')
fig.patch.set_visible(False)
ax.axis('off')
plt.title("E. SST Spatial Gradient, [2006-2026] [$^\circ$C/km]")
#======================================================================
plt.subplot(3,2,6)
#m = Basemap(projection='robin',lon_0=0,resolution='c')
m = Basemap(projection='cyl', llcrnrlon=260, llcrnrlat=7, urcrnrlon=295, urcrnrlat=33)
# m = Basemap(width=6000000,height=4000000,projection='lcc',
#              resolution=None,lat_1=27.,lat_2=30.,lat_0=25.,lon_0=-90.)
# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(LON[0:-1], LAT[0:-1]))
#m.drawcoastlines()
# draw filled contours.
#ax1 = m.contourf(x,y,VAR_DIFF_mask,50,cmap=plt.cm.plasma_r,levels=slevels)
ax1 = m.contourf(x,y,TE_MAG,50,cmap=plt.cm.viridis_r,levels=mlevels)

ax2 = m.contour(x,y,TE_MAG,50,levels=mlevels,colors='k')
plt.clabel(ax2, inline=1, fontsize=10)
cbar = m.colorbar(ax1,location='right',pad="2%")
#m.fillcontinents(color='coral',lake_color='aqua')
cbar.set_label(r'Temperature / km, [$^\circ$C/km]')
fig.patch.set_visible(False)
ax.axis('off')
plt.title("F. SST Spatial Gradient, [2080-2100] [$^\circ$C/km]")


plt.show()



