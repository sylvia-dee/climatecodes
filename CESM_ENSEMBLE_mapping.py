
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

S1=np.load('SALT_ENSMEAN_2006-2080.npy')
S2=np.load('SALT_ENSMEAN_2080-2100.npy')

T1=np.load('SST_ENSMEAN_2006-2080.npy')
T2=np.load('SST_ENSMEAN_2080-2100.npy')

D1=np.load('DIC_ENSMEAN_2006-2080.npy')
D2=np.load('DIC_ENSMEAN_2080-2100.npy')

P1=np.load('PH_ENSMEAN_2006-2080.npy')
P2=np.load('PH_ENSMEAN_2080-2100.npy')

A1=np.load('ALK_ENSMEAN_2006-2080.npy')
A2=np.load('ALK_ENSMEAN_2080-2100.npy')


LAT=np.load('POP2_lats.npy')
LON=np.load('POP2_lon.npy')

LON=LON[0,:]
LAT=LAT[:,0]
#======================================================================
cd /rdf/sd75/sylvia/GULF_CESM/
#======================================================================

#======================================================================
# SAVE AS MEAN FROM 2006-2026, 2080-2100
#======================================================================

PH_1=np.mean(P1[0:240,:,:],axis=0)
PH_2=np.mean(P2,axis=0)

TEMP_1=np.mean(T1[0:240,:,:],axis=0)
TEMP_2=np.mean(T2,axis=0)

VAR_TEMP_1=np.std(np.array(T1[0:240,:,:]),axis=0)
VAR_TEMP_2=np.std(np.array(T2),axis=0)

SALT_1=np.mean(S1[0:240,:,:],axis=0)
SALT_2=np.mean(S2,axis=0)

ALK_1=np.mean(A1[0:240,:,:],axis=0)
ALK_2=np.mean(A2,axis=0)

DIC_1=np.mean(D1[0:240,:,:],axis=0)
DIC_2=np.mean(D2,axis=0)

#longitude=(260.,295.), latitude = (10., 35.))

# lat=np.array(SiO3_80100.getLatitude())[:,0]
# lon=np.array(PO4_0680.getLongitude())[0,:]

#======================================================================
# SAVE FOR MARK

import pandas as pd 
pd.DataFrame(TEMP_1).to_csv("ECESM_GULF_TEMP_2006_2026.csv")
pd.DataFrame(TEMP_2).to_csv("ECESM_GULF_TEMP_2080_2100.csv")

pd.DataFrame(SALT_1).to_csv("ECESM_GULF_SALT_2006_2026.csv")
pd.DataFrame(SALT_2).to_csv("ECESM_GULF_SALT_2080_2100.csv")

pd.DataFrame(PH_1).to_csv("ECESM_GULF_PH_2006_2026.csv")
pd.DataFrame(PH_2).to_csv("ECESM_GULF_PH_2080_2100.csv")

pd.DataFrame(ALK_1).to_csv("ECESM_GULF_ALK_2006_2026_CONCENTRATION.csv")
pd.DataFrame(ALK_2).to_csv("ECESM_GULF_ALK_2080_2100_CONCENTRATION.csv")

pd.DataFrame(DIC_1).to_csv("ECESM_GULF_DIC_2006_2026_CONCENTRATION.csv")
pd.DataFrame(DIC_2).to_csv("ECESM_GULF_DIC_2080_2100_CONCENTRATION.csv")


# pd.DataFrame(LAT).to_csv("ECESM_GULF_LATITUDE.csv")
# pd.DataFrame(LON).to_csv("ECESM_GULF_LONGITUDE.csv")
#======================================================================
#======================================================================
# TEMP
VAR_1_mask=np.ma.masked_outside(VAR_TEMP_1,-10,10)
VAR_2_mask=np.ma.masked_outside(VAR_TEMP_2,-10,10)

#VAR_DIFF=VAR_2_mask-VAR_1_mask
VAR_DIFF=TEMP_2-TEMP_1

VAR_DIFF_mask=np.ma.masked_inside(VAR_DIFF,-0.001,0.001)
#======================================================================

VAR_DIFF=ALK_2-ALK_1
VAR_DIFF_mask=np.ma.masked_inside(VAR_DIFF,-1E-15,1E-15)

MAP_MASK=np.ma.masked_greater(SALT_2,400.)

#MAP=np.ma.set_fill_value(VAR_DIFF,9E9)
#======================================================================

# MAP,lons2 =m.shiftdata(SST_DIFF,Slon, lon_0=0)
from mpl_toolkits.basemap import Basemap, shiftgrid
slevels=np.arange(26,42,0.5)

fig,ax=plt.subplots()
#======================================================================

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
ax1 = m.contourf(x,y,MAP_MASK,50,cmap=plt.cm.BrBG_r,levels=slevels)

ax2 = m.contour(x,y,MAP_MASK,50,levels=slevels,colors='k')
plt.clabel(ax2, inline=1, fontsize=10)
cbar = m.colorbar(ax1,location='bottom',pad="2%")
#m.fillcontinents(color='coral',lake_color='aqua')

fig.patch.set_visible(False)
ax.axis('off')
plt.title("Salinity [2080-2100] [g/kg], CESM RCP8.5 ENSEMBLE")

plt.show()
# #======================================================================
# plt.subplot(1,3,2)
# #m = Basemap(projection='robin',lon_0=0,resolution='c')
# m = Basemap(projection='cyl', llcrnrlon=260, llcrnrlat=10, urcrnrlon=295, urcrnrlat=33)
# # m = Basemap(width=6000000,height=4000000,projection='lcc',
# #              resolution=None,lat_1=27.,lat_2=30.,lat_0=25.,lon_0=-90.)
# # compute the lons and lats to fit the projection
# x, y = m(*np.meshgrid(LON, LAT))
# m.drawcoastlines()
# # draw filled contours.
# ax1 = m.contourf(x,y,VAR_2_mask,50,cmap=plt.cm.Reds,levels=slevels)
# #ax2 = m.contour(x,y,SALT_DIFF,50,levels=slevels,colors='k')
# #plt.clabel(ax2, inline=1, fontsize=10)
# cbar = m.colorbar(ax1,location='bottom',pad="2%")

# fig.patch.set_visible(False)
# ax.axis('off')
# plt.title("SST Standard Dev. 2080-2100")
# #======================================================================
# plt.subplot(1,3,3)
# dlevels=np.arange(-0.2,0.2,0.01)

# #m = Basemap(projection='robin',lon_0=0,resolution='c')
# m = Basemap(projection='cyl', llcrnrlon=260, llcrnrlat=10, urcrnrlon=295, urcrnrlat=33)
# # m = Basemap(width=6000000,height=4000000,projection='lcc',
# #              resolution=None,lat_1=27.,lat_2=30.,lat_0=25.,lon_0=-90.)
# # compute the lons and lats to fit the projection
# x, y = m(*np.meshgrid(LON, LAT))
# m.drawcoastlines()
# # draw filled contours.
# ax1 = m.contourf(x,y,VAR_DIFF,50,cmap=plt.cm.RdBu_r,levels=dlevels)
# #ax2 = m.contour(x,y,SALT_DIFF,50,levels=slevels,colors='k')
# #plt.clabel(ax2, inline=1, fontsize=10)

# fig.patch.set_visible(False)
# ax.axis('off')
# # m.etopo()

# # c,d=m(lons2,lats)
# # cs2 = m.scatter(c,d,c=depth,marker='d',s=75,cmap=plt.cm.plasma,edgecolors='k',linewidth=1.0)


# # add colorbar.
# cbar = m.colorbar(ax1,location='bottom',pad="2%")
# plt.title(r'SST Variance Change [2080-2100] vs.[2006-2026] [C]')

# cbar2 = m.colorbar(cs2,location='right',pad="2%")
# cbar2.set_label(r'Reef Depth [m]')
#plt.title('RCP 8.5 Mean SST 2080-2100 and Gulf Reef System Locations',fontsize=16)
# plt.show()
