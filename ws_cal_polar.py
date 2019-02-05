# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 13:31:20 2018

@author: Jose Luis Rodriguez-Solis
         Oceanografia Fisica
         CICESE

WIND SPEED
Calcula los valores de wind speed basados en :

Archer, C. L., & Caldeira, K. (2008). Historical trends in the jet streams. 
Geophysical Research Letters, 35(8), 1–6. doi.org/10.1029/2008GL033614


"""

from netCDF4 import Dataset

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import numpy as np
#from scipy import stats

import os
import atmos

#import sys

# nivel superior
pi = 400 
# nivel infererior
ps = 100

# donde estan los archivos de Era Interim
path = '/home/luis/Documents/CLASES_CICESE/ERA Interim/datos/mensual01/'
f = os.listdir(path)
f.sort()

ws = []
Ps = []
phi= []
t  = []
ls = []

nc = Dataset (path+f[0])
lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:]
lev = nc.variables['level'][:]

# buscando indices de altura de 400 y 100
zi = sorted(lev).index(pi)
zs = sorted(lev).index(ps)
# Malla para graficar y encontrar latitudes
x,y=np.meshgrid(lon,lat)
# buscando latitudes entre 70N y 15N
yl  = np.where((lat>= 15.0) & (lat<=  60.5)); yl   = yl[0]
ys = yl[0]
yi = yl[-1]

yeari = int(f[0][7:-7])

zl = np.zeros((10,len(lat),len(lon)))
for i in range(0, len(lev[zs:zi+1])):    
    zl[i,:,:] = lev[zs+i]

num = len(f)
for i in range (0,num): # len(f) 300
    
    nc = Dataset (path+f[i])    
    u   = nc.variables['u'][:].squeeze()
    v   = nc.variables['v'][:].squeeze()
    q   = nc.variables['q'][:].squeeze() * 1000.0
    T   = nc.variables['t'][:].squeeze() 
    nc.close()    
    
    # buscando indices de niveles superior e inferior
    us = (u[zs:zi+1,:,:].copy())**2
    del u
    vs = (v[zs:zi+1,:,:].copy())**2
    del v    
    ms = np.sqrt(us + vs)
    qs = q[zs:zi+1,:,:].copy() 
    del q 
    ts = T[zs:zi+1,:,:].copy()
    del T

    # calculando rho
    rho = atmos.calculate('rho',p = zl*100.0, Tv = ts)
    
    # calculando el wind speed
    ws.append ( (np.nansum( rho * ms, axis=0)) / (np.nansum(rho, axis=0) ) )
    # calculando Presion
    Ps.append ( (np.nansum( rho * ms * zl, axis=0 )) / (np.nansum(rho * ms, axis=0) ) )
    # calculando latitud
    ls_top = np.nansum ( rho * ms, axis=0) * y
    ls_top = np.nansum ( ls_top[ys:yi,:], axis = 0 )
    
    ls_inf = np.nansum ( rho * ms, axis=0)
    ls_inf = np.nansum ( ls_inf[ys:yi,:], axis = 0 )
    
    ls.append( ls_top/ls_inf )
    
    t.append(f[i][7:-3])
    print f[i][7:-3]

yearf = int(f[i][7:-7])

print ('finaliza ciclo de calculos')
print ('')
print ('comienza ciclo para medias')
#%% inicia promedios 

ws = np.asarray(ws)
wsm = np.nanmean(ws,axis = 0)

wn = np.nanmean (ws[:,ys:yi,:], axis = 1)
wn = np.nanmean (wn, axis = 1)

Ps = np.asarray(Ps)
psm = np.nanmean(Ps,axis = 0)

pn = np.nanmean (Ps[:,ys:yi,:], axis = 1)
pn = np.nanmean (pn, axis = 1)

ls = np.asarray(ls)
ls = np.nanmean(ls,axis=1)
#%% VIENTO
# promedios para verano
print ('viento')
wn_ver = []
wnsver = []
k = 0
for j in range(5,len(pn),12):
    print t[j:j+3]
    wn_ver.append(np.nanmean (wn[j:j+3]))
    wnsver.append(np.nanmean (ws[j:j+3,:,:],axis=0 ) )
wnsver = np.asarray( wnsver )
# promedios para invierno
wn_inv = []
wnsinv = []
k = 0
for j in range(11,len(wn)-1,12):
    print t[j:j+3]
    wn_inv.append(np.nanmean (wn[j:j+3]))    
    wnsinv.append(np.nanmean (ws[j:j+3,:,:],axis=0 ) )
wnsinv = np.asarray( wnsinv )
# promedio anual    
wni = []
k = 0
for i in range(12,num+1,12):
    print i
    wni.append( np.nanmean(ws[k:i,ys:yi+1,:] ) )
    k = k + 12
wni = np.asarray(wni)
    
#%% PRESION
# promedios para verano
print ('presion')
pn_ver = []
psver  = []
k = 0
for j in range(5,len(pn),12):
    print t[j:j+3]
    pn_ver.append(np.nanmean (pn[j:j+3]))
    psver.append ( np.nanmean(Ps[j:j+3,:,:],axis = 0) )
psver = np.asarray(psver)
# promedios para invierno
pn_inv = []
psinv  = []
k = 0
for j in range(11,len(pn)-1,12):
    print t[j:j+3]
    pn_inv.append(np.nanmean (pn[j:j+3]))    
    psinv.append ( np.nanmean(Ps[j:j+3,:,:],axis = 0) )
psinv = np.asarray(psinv)
# promedio anual    
pni = []
k = 0
for i in range(12,num+1,12):
    print i
    pni.append( np.nanmean(Ps[k:i,ys:yi+1,:] ) )
    k = k + 12
pni = np.asarray(pni)    
    
#%% LATITUD
print ('Latitud')
# promedios para latitud anual
ln = []
k = 0
for j in range(0,len(pn),12):
    print t[j:j+12]
    ln.append(np.nanmean (ls[j:j+12]))    
ln = np.asarray(ln)
# promedios para verano
ln_ver = []
k = 0
for j in range(5,len(pn),12):
    print t[j:j+3]
    ln_ver.append(np.nanmean (ls[j:j+3]))
ln_ver = np.asarray( ln_ver )
# promedios para invierno
ln_inv = []
k = 0
for j in range(11,len(wn)-1,12):
    print t[j:j+3]
    ln_inv.append(np.nanmean (ls[j:j+3]))    
ln_inv = np.asarray( ln_inv )
#%%
ya = -90           # latitud inicial
yf =  90            # latitud final
xa = 0           # longitud inicial
xf = 360           # longitud final

#%% promedio anual viento
print ('comienza a graficar')


#%%

clev = np.arange(10,61,5)

fig, ax = plt.subplots(num=1, figsize=(12,5))

plt.subplot(131)
m = Basemap(projection='npstere',boundinglat=0,lon_0=260,resolution='l')
#m = Basemap(lon_0=-105,lat_0=90,projection='ortho')
xpo, ypo =m(x,y)
m.drawcoastlines(linewidth=0.25)
m.drawcountries(linewidth=0.25)
m.drawparallels(np.arange(-80.,81.,20.),linewidth=0.25,labels=[1,0,1,1],fontsize=8)
m.drawmeridians(np.arange(-180.,181.,20.),linewidth=0.25,labels=[1,0,1,1],fontsize=8) 
cs = m.contourf(xpo,ypo,wsm,levels=clev,cmap=plt.cm.Greens,extend='both',round=True)
plt.text (0.1,0.1,u'Anual',transform=ax.transAxes,fontsize=10,color='k')

plt.subplot(132)
m = Basemap(projection='npstere',boundinglat=0,lon_0=260,resolution='l')
xpo, ypo =m(x,y)
m.drawcoastlines(linewidth=0.25)
m.drawcountries(linewidth=0.25)
m.drawparallels(np.arange(-80.,91.,20.),linewidth=0.25,labels=[0,0,1,1],fontsize=8)
m.drawmeridians(np.arange(-180.,181.,20.),linewidth=0.25,labels=[0,0,1,1],fontsize=8) 
cs = m.contourf(xpo,ypo,np.nanmean(wnsinv,axis=0),levels=clev,cmap=plt.cm.Greens,extend='both',round=True)
plt.title ('Promedio de velocidad del viento (Jet Stream)\nDatos:Era Interim 1981-2017\n',fontsize=10,loc = 'center')
plt.text (0.43,0.1,u'Invierno',transform=ax.transAxes,fontsize=10,color='k')

plt.subplot(133)
m = Basemap(projection='npstere',boundinglat=0,lon_0=260,resolution='l')
xpo, ypo =m(x,y)
m.drawcoastlines(linewidth=0.25)
m.drawcountries(linewidth=0.25)
m.drawparallels(np.arange(-80.,81.,20.),linewidth=0.25,labels=[0,1,1,1],fontsize=8)
m.drawmeridians(np.arange(-180.,181.,20.),linewidth=0.25,labels=[0,1,1,1],fontsize=8) 
cs = m.contourf(xpo,ypo,np.nanmean(wnsver,axis=0),levels=clev,cmap=plt.cm.Greens,extend='both',round=True)
plt.text (0.75,0.1,u'Verano',transform=ax.transAxes,fontsize=10,color='k')

fig.subplots_adjust(right=0.82)
#[left, bottom, width, height] 
cbar_ax = fig.add_axes([0.85, 0.25, 0.009, 0.5])
cbar = fig.colorbar(cs, cax=cbar_ax)
cbar.set_label(u'm/s',fontsize=8)
plt.savefig('PolarWS.png',bbox_inches='tight',dpi=150)


#%%
clev = np.arange(230,301,5)

fig, ax = plt.subplots(num=2, figsize=(12,5))

plt.subplot(131)
m = Basemap(projection='npstere',boundinglat=0,lon_0=260,resolution='l')
xpo, ypo =m(x,y)
m.drawcoastlines(linewidth=0.25)
m.drawcountries(linewidth=0.25)
m.drawparallels(np.arange(-80.,81.,20.),linewidth=0.25,labels=[1,0,1,1],fontsize=8)
m.drawmeridians(np.arange(-180.,181.,20.),linewidth=0.25,labels=[1,0,1,1],fontsize=8) 
cs = m.contourf(xpo,ypo,psm,levels=clev,cmap=plt.cm.bone_r,extend='both',round=True)
plt.text (0.1,0.1,u'Anual',transform=ax.transAxes,fontsize=10,color='k')

plt.subplot(132)
m = Basemap(projection='npstere',boundinglat=0,lon_0=260,resolution='l')
xpo, ypo =m(x,y)
m.drawcoastlines(linewidth=0.25)
m.drawcountries(linewidth=0.25)
m.drawparallels(np.arange(-80.,91.,20.),linewidth=0.25,labels=[0,0,1,1],fontsize=8)
m.drawmeridians(np.arange(-180.,181.,20.),linewidth=0.25,labels=[0,0,1,1],fontsize=8) 
cs = m.contourf(xpo,ypo,np.nanmean(psinv,axis=0),levels=clev,cmap=plt.cm.bone_r,extend='both',round=True)
plt.title (u'Promedio de Presión (Jet Stream)\nDatos:Era Interim 1981-2017\n',fontsize=10,loc = 'center')
plt.text (0.43,0.1,u'Invierno',transform=ax.transAxes,fontsize=10,color='k')

plt.subplot(133)
m = Basemap(projection='npstere',boundinglat=0,lon_0=260,resolution='l')
xpo, ypo =m(x,y)
m.drawcoastlines(linewidth=0.25)
m.drawcountries(linewidth=0.25)
m.drawparallels(np.arange(-80.,81.,20.),linewidth=0.25,labels=[0,1,1,1],fontsize=8)
m.drawmeridians(np.arange(-180.,181.,20.),linewidth=0.25,labels=[0,1,1,1],fontsize=8) 
cs = m.contourf(xpo,ypo,np.nanmean(psver,axis=0),levels=clev,cmap=plt.cm.bone_r,extend='both',round=True)
plt.text (0.75,0.1,u'Verano',transform=ax.transAxes,fontsize=10,color='k')
#cbar = m.colorbar(cs,location='right',ticks=clev,size="05%",pad="15%",extend='both',drawedges=True)
#cbar.set_label(u'm/s',x=1.05,fontsize=8)

fig.subplots_adjust(right=0.82)
#[left, bottom, width, height] 
cbar_ax = fig.add_axes([0.85, 0.25, 0.009, 0.5])
cbar = fig.colorbar(cs, cax=cbar_ax)
cbar.set_label(u'hPa',fontsize=8)
plt.savefig('PolarPS.png',bbox_inches='tight',dpi=150)


