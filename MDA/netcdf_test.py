#!/usr/bin/python3

from scipy.io import netcdf_file
#name = '05_Prod4_tgsn166.nc'
name = "/home/hgibbs/research/moore/vpu/04_Prod_v21l.nc"
f=netcdf_file(name,'r')
print(f.application)
coords = f.variables['coordinates']
print(coords.units)
print(coords.shape)

data0=coords[0].copy()
print(coords[-1])
f.close()

