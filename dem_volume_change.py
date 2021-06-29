#!/usr/bin/python
#
# calculate glacier volume change from 2 DEMs of different years
# lauren vargo, august 2020
# uses conda environment, conda activate /Users/home/vargola/geo4520 before using
# helpful: https://pybob.readthedocs.io/en/stable/modules/index.html

from pybob.GeoImg import GeoImg
import pybob.ddem_tools as ddem
import pybob.coreg_tools as ct
import pybob.image_tools as imtool
import pybob.coreg_tools as coreg
from pybob.plot_tools import set_pretty_fonts
import matplotlib.pyplot as plt
import numpy as np
import argparse
from osgeo import gdal
import geopandas as gpd
from PIL import Image
import os

### parser arguments
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('dem_1', help='dem of earlier year')
parser.add_argument('dem_2', help='dem of later year')
args = parser.parse_args()
coregFlag = True   # True includes coregistration
glacFlag = 'Brewster' # 'Rolleston or 'Brewster'
if glacFlag == 'Rolleston':
    odir = '/Volumes/arc_03/vargola/eoss_images/ROLLESTON/'
    bedrx_poly_mask = odir+'bedrx_mask3.shp'
else: 
    odir = '/Volumes/arc_03/vargola/eoss_images/BREWSTER/'
    bedrx_poly_mask = odir+'bedrx_poly.shp'

# read in DEMs and project to same extents
dem1 = GeoImg(args.dem_1) 
dem2 = GeoImg(args.dem_2)
x_res = dem1.dx  # had to add this for calculating dv at the end, mask seems to distort 
y_res = dem1.dy

#coregister
if coregFlag:
    coreg.dem_coregistration(dem1,dem2,glaciermask=None,landmask=bedrx_poly_mask,outdir=odir,pts=False)
    dem2 = GeoImg(args.dem_2[0:-4] + '_adj.tif')

# dh_dem
dem2_reproj = dem2.reproject(dem1) # reproject master DEM to slave DEM extent, cell size
dh_dem = dem2.copy(new_raster=dem2_reproj.img - dem1.img)
dh_bedrx = dem2.copy(new_raster=dem2_reproj.img - dem1.img)

# remove outliers
dem_rem = np.absolute(np.nanmean(dh_dem.img))+(np.nanstd(dh_dem.img)*5)
np.seterr(all='ignore')
dh_dem.img[dh_dem.img > dem_rem] = np.nan
dh_dem.img[dh_dem.img < (dem_rem*-1)] = np.nan
dh_bedrx.img[dh_bedrx.img > dem_rem] = np.nan
dh_bedrx.img[dh_bedrx.img < (dem_rem*-1)] = np.nan

# test coregistered bedrx error
bedrx_ma = imtool.create_mask_from_shapefile(dem1,bedrx_poly_mask,buffer=None)
bedrx_mask = np.invert(bedrx_ma) # inverting for mask
dh_bedrx.img = np.ma.masked_where(bedrx_mask, dh_bedrx.img) # masking
fig0 = dh_bedrx.display(fig=plt.figure(figsize=(8,8)), cmap='bwr_r', vmin=-10, vmax=10)
fig0.gca().set_xlabel('utm easting (m)')
fig0.gca().set_ylabel('utm northing (m)')
cb = plt.colorbar(); cb.set_label('elevation difference (m)')
#ts = np.ma.compressed(dh_bedrx.img)
#plt.hist(ts,bins=40,range=[-7.5,7.5])
#plt.show()  # to plot histogram

# mask DEMs
s1 = gpd.read_file(odir+'2016_outline1.shp')
s2 = gpd.read_file(odir+'2021_outline.shp')
s1arr_a = imtool.create_mask_from_shapefile(dem1,odir+'2016_outline1.shp',buffer=None)
s1arr = np.invert(s1arr_a)
dh_dem.img = np.ma.masked_where(s1arr, dh_dem.img) # masking

# plot
fig1 = dh_dem.display(fig=plt.figure(figsize=(8,8)), cmap='bwr_r', vmin=-10, vmax=10)
fig1.gca().set_xlabel('utm easting (m)')
fig1.gca().set_ylabel('utm northing (m)')
cb = plt.colorbar(); cb.set_label('elevation difference (m)')

# calculate volume change
dv = np.nansum(dh_dem.img) * np.abs(x_res) * np.abs(y_res)
s1_area = np.nansum(s1.area)
s2_area = np.nansum(s2.area)
smean = (s2_area+s1_area)/2
density = 0.85  # 0.85 p/m 0.06 in Huss density paper, 0.92 google 
B = (dv/smean) * (density/1)
print(B)

# uncertainties, all in %
sig_geod_obs = (np.nanstd(dh_bedrx.img)/np.absolute(np.nanmean(dh_dem.img)))*100  # DEM uncertainty
sig_calib = 0 
sig_density = (60/850)*100 # 850 p/m 60 kg m-2
sig_geod = (sig_geod_obs**2 + sig_calib**2 + sig_density**2)**0.5
sig_glac = 0 
sig_area = 5  # huss hock uses 5%, we could estimate our own, should be lower
sig = (sig_glac**2 + sig_geod**2 + sig_area**2)**0.5
sig = B * (sig*0.01)  # need to convert back from percent
print(sig)
