

import numpy as np
import numpy.ma as ma

import string
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import sys, os, cdtime, cdutil, cdms2, MV2, time, datetime
from math import exp
import scipy.io as sio
#======================================================================


cd /rdf/sd75/sylvia/GULF_CESM/
#======================================================================


T1=np.load('SST_ENSMEAN_2006-2080.npy')
T2=np.load('SST_ENSMEAN_2080-2100.npy')

LAT=np.load('POP2_lats.npy')
LON=np.load('POP2_lon.npy')

LON=LON[0,:]
LAT=LAT[:,0]

cd /rdf/sd75/sylvia/GOM/

import numpy as np
import csv
from numpy import genfromtxt

data = genfromtxt('Corals_2.csv',
                delimiter = ',',
                names=True,
                dtype=None
                #mask module may be used when better understood
                #,usemask=True
                )
depth=data['depth']
lats=data['lat']
lons=data['lon']


#======================================================================
# We computed the cumulative SST months above the mean of maximum monthly 
# climatology in each Gulf of Mexico grid cell 

# DHM threshold exceeding 1$^\circ$C in a given month will lead to bleaching; 
# consider temperatures exceeding 2$^\circ$C above monthly climatologies are considered
# a higher threshold, corresponding to a degree heating week (DHW) exceeding 8 weeks 
# of high heating. 
#======================================================================

# COMPUTE MAXIMUM MONTHLY CLIMATOLOGY (warmest month in the year)

TEMP_1=T1
TEMP_2=T2

seasons1=np.reshape(TEMP_1,[75,12,384, 320])
seasons1=np.ma.masked_greater(seasons1, 500)
climat=np.mean(seasons1,axis=0)

hotmonth=np.max(climat,axis=0)

# compute number of months above DHM threshold of 1, 2. 

DHW_1 =np.zeros((T1.shape))
DHW_2 =np.zeros((T1.shape))

for u in range(240):
	for i in range(len(LAT)):
		for j in range(len(LON)):
			if TEMP_1[u,i,j] > (hotmonth[i,j]+1):
				DHW_1[u,i,j] = 1.0
			else:
				DHW_1[u,i,j] = 0.0

#======================================================================

for u in range(240):
	for i in range(len(LAT)):
		for j in range(len(LON)):
			if TEMP_2[u,i,j] > (hotmonth[i,j]+1):
				DHW_2[u,i,j] = 1.0
			else:
				DHW_2[u,i,j] = 0.0

#======================================================================
cd /rdf/sd75/sylvia/GULF_CESM/


np.save('DHW1.npy',DHW_1)
np.save('DHW2.npy',DHW_2)


DHW_1=np.load('DHW1.npy')
DHW_2=np.load('DHW2.npy')


DHW_TOT_1=np.sum(DHW_1,axis=0)
DHW_TOT_2=np.sum(DHW_2,axis=0)

DIFF=DHW_TOT_2-DHW_TOT_1

D1_mask=np.ma.masked_inside(DHW_TOT_1,-0.001,0.001)
D2_mask=np.ma.masked_inside(DHW_TOT_2,-0.001,0.001)
DIFF_mask=np.ma.masked_inside(DIFF,-0.001,0.001)
#=====================================================================

#MAP=np.ma.set_fill_value(VAR_DIFF,9E9)
#======================================================================

# MAP,lons2 =m.shiftdata(SST_DIFF,Slon, lon_0=0)
from mpl_toolkits.basemap import Basemap, shiftgrid
slevels=np.arange(40,100,5)

fig,ax=plt.subplots()
#======================================================================

plt.subplot(1,3,1)
#m = Basemap(projection='robin',lon_0=0,resolution='c')
m = Basemap(projection='cyl', llcrnrlon=260, llcrnrlat=7, urcrnrlon=295, urcrnrlat=33)
# m = Basemap(width=6000000,height=4000000,projection='lcc',
#              resolution=None,lat_1=27.,lat_2=30.,lat_0=25.,lon_0=-90.)
# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(LON, LAT))
#m.drawcoastlines()
# draw filled contours.
#ax1 = m.contourf(x,y,VAR_DIFF_mask,50,cmap=plt.cm.plasma_r,levels=slevels)
ax1 = m.contourf(x,y,DHW_TOT_1,50,cmap=plt.cm.plasma_r,levels=slevels)
ax2 = m.contour(x,y,DHW_TOT_1,50,levels=slevels,colors='k')

#plt.clabel(ax2, inline=1, fontsize=10)
cbar = m.colorbar(ax1,location='bottom',pad="2%")
cbar.set_label(r'DEGREE HEATING MONTHS (DHM)',fontsize=14)

fig.patch.set_visible(False)
ax.axis('off')
plt.title("DHM 2006-2026",fontsize=16)
#======================================================================

plt.subplot(1,3,2)
#m = Basemap(projection='robin',lon_0=0,resolution='c')
m = Basemap(projection='cyl', llcrnrlon=260, llcrnrlat=7, urcrnrlon=295, urcrnrlat=33)
# m = Basemap(width=6000000,height=4000000,projection='lcc',
#              resolution=None,lat_1=27.,lat_2=30.,lat_0=25.,lon_0=-90.)
# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(LON, LAT))
#m.drawcoastlines()
# draw filled contours.
#ax1 = m.contourf(x,y,VAR_DIFF_mask,50,cmap=plt.cm.plasma_r,levels=slevels)
ax1 = m.contourf(x,y,D2_mask,50,cmap=plt.cm.plasma_r,levels=slevels)
ax2 = m.contour(x,y,D2_mask,50,levels=slevels,colors='k')
#ax2 = m.contour(x,y,MAP_MASK,50,levels=slevels,colors='k')
#plt.clabel(ax2, inline=1, fontsize=10)
cbar = m.colorbar(ax1,location='bottom',pad="2%")
cbar.set_label(r'DEGREE HEATING MONTHS (DHM)',fontsize=14)

fig.patch.set_visible(False)
ax.axis('off')
plt.title("DHM 2080-2100",fontsize=16)

#======================================================================
fig,ax=plt.subplots()
plt.subplot(1,1,1)
#m = Basemap(projection='robin',lon_0=0,resolution='c')
m = Basemap(projection='cyl', llcrnrlon=260, llcrnrlat=7, urcrnrlon=295, urcrnrlat=33)
# m = Basemap(width=6000000,height=4000000,projection='lcc',
#              resolution=None,lat_1=27.,lat_2=30.,lat_0=25.,lon_0=-90.)
# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(LON, LAT))
#m.drawcoastlines()
# draw filled contours.
#ax1 = m.contourf(x,y,VAR_DIFF_mask,50,cmap=plt.cm.plasma_r,levels=slevels)
ax1 = m.contourf(x,y,DIFF_mask,50,cmap=plt.cm.plasma_r,levels=slevels)
ax2 = m.contour(x,y,DIFF_mask,50,levels=slevels,colors='k')
#ax2 = m.contour(x,y,MAP_MASK,50,levels=slevels,colors='k')
plt.clabel(ax2, inline=1, fontsize=12,fmt='%1.0f')
cbar = m.colorbar(ax1,location='bottom',pad="2%")
cbar.set_label(r'# of Degree Heating Months (DHM), difference in total month count',fontsize=14)
#m.fillcontinents(color='coral',lake_color='aqua')

fig.patch.set_visible(False)
ax.axis('off')
plt.title("DHM Late 21st Cen. Difference [2080-2100] minus [2006-2026]",fontsize=16)

plt.show()

