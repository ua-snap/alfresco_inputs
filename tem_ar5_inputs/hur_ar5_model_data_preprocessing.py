# now lets try xray
import pandas as pd
import numpy as np
import os, sys, re, xray
from rasterio import Affine as A
from rasterio.warp import reproject, RESAMPLING
from osgeo import gdal
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid

# some setup pathing
input_dir = '~/Documents/hur'
os.chdir( input_dir )

# the level of the atmosphere we want to use
atmos_level = 11

# open multiple datasets as a single file
xds = xray.open_mfdataset( 'hur_Amon_GFDL-CM3_historical_r1i1p1_*.nc' )
xds_hur = xds.hur.loc['1900-01-01':'2005-12-12'] # slice the dataset using the time variable in xray object
hur_lev = xds_hur[ :, atmos_level, ... ]

# calculate climatology and anomalies
climatology = hur_lev.loc[ '1961-01-01':'1990-12-31' ].groupby( 'time.month' ).mean( 'time' )
anomalies = hur_lev.groupby( 'time.month' ) - climatology

# # # REPROJECT AND CROP EXTENT
# what do we need to do to properly resample the data
time_len, rows, cols = hur_lev.shape
# NOTE: geotransform = [left, res, 0.0, top, 0.0, res]
height = rows
width = cols
crs = 'epsg:4326'
affine = A( *[np.diff( xds.lon )[ 0 ], 0.0, -180.0, 0.0, -np.diff( xds.lat )[ 0 ], 90.0] )
count = time_len
resolution = ( np.diff( xds.lat )[ 0 ], np.diff( xds.lon )[ 0 ] )

# use basemap's shiftgrid to shift to the needed -180to180 longitude ordering
# dat3, lons = shiftgrid( 180., hur_lev.values[0,...][:], hur_lev.lon.data, start=False )#

# lets try to do this across all bands
# test = np.apply_over_axes( lambda x: shiftgrid( 180., x, a=hur_lev.lon.data, start=False ), hur_lev.values, axes=(2,3) )
# this works for shifting the grid across multiple bands.
dat, lons = shiftgrid( 180., anomalies[:], anomalies.lon.data, start=False )#

# metadata for input?
meta = { 'affine':affine,
		'height':rows,
		'width':cols,
		'crs':crs,
		'driver':'GTiff',
		'dtype':np.float32,
		'count':1 }

meta.update( transform=meta['affine'].to_gdal() )

# test write
output_filename = '/home/UA/malindgren/Documents/test_out_hur_orient.tif'
with rasterio.open( output_filename, 'w', **meta ) as out:
	out.write_band( 1, hur_lev[1, ...].astype( np.float32 ) )

# test write
output_filename = '/home/UA/malindgren/Documents/test_out_rst_c.tif'
with rasterio.open( output_filename, 'w', **meta ) as out:
	out.write_band( 1, dat.astype( np.float32 ) )

# setup output and reproject
rst = rasterio.open( '/home/UA/malindgren/Documents/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif' )
dst = np.empty_like( rst.read( 1 ) )

reproject( dat,
			dst,
			src_transform=affine.to_gdal(),
			src_crs={'init':'epsg:4326'},
			dst_transform=rst.affine.to_gdal(),
			dst_crs={'init':'epsg:3338'},
			resampling=RESAMPLING.cubic_spline )

# write it out
output_filename = '/home/UA/malindgren/Documents/test_out_rst5.tif'
meta = rst.meta
meta.update( compress='lzw' ) #, crs={'init':'epsg:3338'}
with rasterio.open( output_filename, 'w', **meta ) as out:
	out.write_band( 1, dst.astype( np.float32 ) )

# plot it with mpl
plt.imshow( dst )
plt.savefig( '/home/UA/malindgren/Documents/test_out5.png' )
plt.close()
