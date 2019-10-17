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

cd /rdf/sd75/sylvia/LMR/

S=sio.loadmat('mnsst_diff_2080-2100_minus_1920-2005_new.mat')
latS=S['tlat']
lonS=S['tlon']
sst_20th=S['sstmn20th'] #1920 to 2005
sst_21st=S['sstmn21st'] # 2080 to 2100

Slon=lonS[0,:]
Slat=latS[:,0]

#stoney corals, gorgonian corals and sponges
# #======================================================================
cd /rdf/sd75/sylvia/GOM/

import numpy as np
import csv
from numpy import genfromtxt

#data = genfromtxt('Corals_2.csv',
data =genfromtxt('Corals_REV1.csv',
                delimiter = ',',
                names=True,
                dtype=None
                #mask module may be used when better understood
                #,usemask=True
                )
depth=-data['D']
lats=data['LAT']
lons=data['LON']
species=data['SP']

lons2=lons+360.

blats=[]
blons=[]
bd=[]

slats=[]
slons=[]
sd=[]

glats=[]
glons=[]
gd=[]

stlats=[]
stlons=[]
std=[]

plats=[]
plons=[]
pd=[]

llats=[]
llons=[]
ld=[]

olats=[]
olons=[]
od=[]

for i in range(len(lons)):
    if species[i]=="BC":
        blats.append(lats[i])
        blons.append(lons[i])
        bd.append(depth[i])

for i in range(len(lons)):
    if species[i]=="SPONGE":
        slats.append(lats[i])
        slons.append(lons[i])
        sd.append(depth[i])

for i in range(len(lons)):
    if species[i]=="GORGO":
        glats.append(lats[i])
        glons.append(lons[i])
        gd.append(depth[i])

for i in range(len(lons)):
    if species[i]=="STONY":
        stlats.append(lats[i])
        stlons.append(lons[i])
        std.append(depth[i])

for i in range(len(lons)):
    if species[i]=="PEN":
        plats.append(lats[i])
        plons.append(lons[i])
        pd.append(depth[i])

for i in range(len(lons)):
    if species[i]=="LACE":
        llats.append(lats[i])
        llons.append(lons[i])
        ld.append(depth[i])

for i in range(len(lons)):
    if species[i]=="SOFT":
        olats.append(lats[i])
        olons.append(lons[i])
        od.append(depth[i])

# for i in range(len(species)):
# 	if species[i]=="BC":
# 		markers[i]='d'
# 	elif species[i]=="SPONGE":
# 		markers[i]='>'
# 	elif species[i]=="GORGO":
# 		markers[i]='s'
# 	elif species[i]=="STONY":
# 		markers[i]='o'
# 	elif species[i]=="PEN":
# 		markers[i]='2'
# 	elif species[i]=="LACE":
# 		markers[i]='x'
# 	elif species[i]=="SOFT":
# 		markers[i]='p'

#======================================================================

# SST_DIFF=sst_21st-sst_20th

slevels=np.arange(25,35,0.5)
smax=35
smin=25

# MAP,lons2 = shiftgrid(350.,SST_DIFF,Slon,start=False)

# MAP,lons2 =m.shiftdata(SST_DIFF,Slon, lon_0=0)

#======================================================================

#m = Basemap(projection='robin',lon_0=0,resolution='c')
#m = Basemap(projection='cyl', llcrnrlon=260, llcrnrlat=10, urcrnrlon=295, urcrnrlat=35)

# # draw costlines and coutries
#m.drawcoastlines(linewidth=1.5)
#m.drawcountries(linewidth=1.5)

from mpl_toolkits.basemap import Basemap, shiftgrid
#======================================================================

plt.subplot(1,2,1)
m = Basemap(width=5000000,height=3000000,projection='lcc',
            resolution=None,lat_1=25.,lat_2=26.,lat_0=22.,lon_0=-85.)

#x, y = m(*np.meshgrid(Slon, Slat))

m.bluemarble()
# ax1 = m.contourf(x,y,sst_21st,50,cmap=plt.cm.plasma_r,levels=slevels,vmin=smin, vmax=smax)
# ax2 = m.contour(x,y,sst_21st,20,levels=slevels,colors='k')
# plt.clabel(ax2, inline=1, fontsize=10,fmt='%2.1f')


#======================================================================
a,b=m(blons,blats)
c,d=m(slons,slats)
o,p=m(glons,glats)
e,f=m(stlons,stlats)
g,h=m(plons,plats)
k,l=m(llons,llats)
u,w=m(olons,olats)
#======================================================================

#cs2 = m.scatter(a,b,c=bd,marker='d',s=75, cmap=plt.cm.YlGnBu_r,edgecolors='k',linewidth=1.0,label='black coral')
#cs3 = m.scatter(c,d,c=sd,marker='>',s=75, cmap=plt.cm.YlGnBu_r,edgecolors='k',linewidth=1.0,label='sponge')
#cs4 = m.scatter(o,p,c=gd,marker='s',s=75, cmap=plt.cm.YlGnBu_r,edgecolors='k',linewidth=1.0,label='gorgonian')
cs5 = m.scatter(e,f,c=std,marker='o',s=75, cmap=plt.cm.YlGnBu_r,edgecolors='k',linewidth=1.0,label='Stony Coral')
#cs6 = m.scatter(g,h,c=pd,marker='2',s=75, cmap=plt.cm.YlGnBu_r,edgecolors='k',linewidth=1.0,label='sea pen')
#cs7 = m.scatter(k,l,c=ld,marker='x',s=75, cmap=plt.cm.YlGnBu_r,edgecolors='k',linewidth=1.0,label='lace coral')
#cs8 = m.scatter(u,w,c=od,marker='p',s=75, cmap=plt.cm.YlGnBu_r,edgecolors='k',linewidth=1.0,label='soft coral')

#cbar_s = m.colorbar(cs2,location='right',pad="2%")
#plt.legend(loc='upper center', bbox_to_anchor=(0.0, -0.3), shadow=True, ncol=1)
plt.legend(loc='upper center', ncol=1, fontsize=12,bbox_to_anchor=(0.5, -0.05))

#======================================================================
# add colorbar.
# cbar = m.colorbar(ax1,location='bottom',pad="2%")
# cbar.set_label(r'SST Mean [C]')

cbar2 = m.colorbar(cs2,location='right',pad="2%")
cbar2.set_label(r'Reef Depth [m]')
#plt.title('RCP 8.5 Mean SST 2080-2100 and Gulf Reef System Locations',fontsize=16)
plt.title('A. Gulf Stony Coral locations',fontsize=16)
#======================================================================

plt.subplot(1,2,2)
m = Basemap(width=5000000,height=3000000,projection='lcc',
            resolution=None,lat_1=25.,lat_2=26.,lat_0=22.,lon_0=-85.)

#x, y = m(*np.meshgrid(Slon, Slat))

m.bluemarble()
# ax1 = m.contourf(x,y,sst_21st,50,cmap=plt.cm.plasma_r,levels=slevels,vmin=smin, vmax=smax)
# ax2 = m.contour(x,y,sst_21st,20,levels=slevels,colors='k')
# plt.clabel(ax2, inline=1, fontsize=10,fmt='%2.1f')

#======================================================================
a,b=m(blons,blats)
c,d=m(slons,slats)
o,p=m(glons,glats)
e,f=m(stlons,stlats)
g,h=m(plons,plats)
k,l=m(llons,llats)
u,w=m(olons,olats)
#======================================================================

cs2 = m.scatter(a,b,c=bd,marker='d',s=75, cmap=plt.cm.YlGnBu_r,edgecolors='k',linewidth=1.0,label='black coral')
cs3 = m.scatter(c,d,c=sd,marker='>',s=75, cmap=plt.cm.YlGnBu_r,edgecolors='k',linewidth=1.0,label='sponge')
cs4 = m.scatter(o,p,c=gd,marker='s',s=75, cmap=plt.cm.YlGnBu_r,edgecolors='k',linewidth=1.0,label='gorgonian')
#cs5 = m.scatter(e,f,c=std,marker='o',s=75, cmap=plt.cm.YlGnBu_r,edgecolors='k',linewidth=1.0,label='stony coral')
cs6 = m.scatter(g,h,c=pld,marker='2',s=75, cmap=plt.cm.YlGnBu_r,edgecolors='k',linewidth=1.0,label='sea pen')
cs7 = m.scatter(k,l,c=ld,marker='x',s=75, cmap=plt.cm.YlGnBu_r,edgecolors='k',linewidth=1.0,label='lace coral')
cs8 = m.scatter(u,w,c=od,marker='p',s=75, cmap=plt.cm.YlGnBu_r,edgecolors='k',linewidth=1.0,label='soft coral')

#cbar_s = m.colorbar(cs2,location='right',pad="2%")
#plt.legend(loc='upper center', bbox_to_anchor=(0.0, -0.3), shadow=True, ncol=1)
plt.legend((cs2,cs3,cs4,cs6,cs7,cs8),('Black Coral','Sponge','Gorgonian','Sea Pen','Lace Coral','Soft Corals'),numpoints=1, loc='upper center', ncol=3, fontsize=12,bbox_to_anchor=(0.5, -0.05))
#======================================================================
# add colorbar.
# cbar = m.colorbar(ax1,location='bottom',pad="2%")
# cbar.set_label(r'SST Mean [C]')

cbar2 = m.colorbar(cs2,location='right',pad="2%")
cbar2.set_label(r'Reef Depth [m]')
#plt.title('RCP 8.5 Mean SST 2080-2100 and Gulf Reef System Locations',fontsize=16)
plt.title('B. Gulf Reefs: all other species',fontsize=16)



plt.show()