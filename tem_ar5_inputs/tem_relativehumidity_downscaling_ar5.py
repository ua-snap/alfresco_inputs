# # # 
# some new code for processing netcdf cmip5 data in python for downscaling
# # # 

class RelativeHumidityDS( object ):
	'''
	class to store information and the data for a specific 
	experiment (historical, rcp26, etc...) stacked into a single array 
	from a potential set of netcdf files.

	takes as input:

	input_dir = string path to a directory storing a single set of 
		model / variable / experiment / realization physics / ... 
	var = string abbreviated name for the variable i.e. "hur", "tas", "rsds"
	level = integer value used to get a specific level from the data set 

	'''
	def __init__( self, input_dir, var, level, begin_subset_date, end_subset_date, climatology_begin_date, climatology_end_date, **kwargs ):
		# filelisting
		filelist = glob.glob( os.path.join( input_dir, '*.nc' ) )
		filelist.sort()
		self.filelist = filelist
		self.var = var
		self.level = level
		self.begin_subset_date = begin_subset_date
		self.end_subset_date = end_subset_date
		dates, data = self.get_subset_data()
		self.data = data
		self.dates = dates
		self.shape = self.data.shape
		self.climatology_begin_date = climatology_begin_date
		self.climatology_end_date = climatology_end_date
	def get_metadata( self ):
		None
	def get_data( self, **kwargs ):
		from scipy.io import netcdf_file
		lev, = np.where( netcdf_file( self.filelist[0] ).variables[ 'plev' ].data == self.level )
		lev = lev.tolist()[ 0 ] 
		return np.vstack([ netcdf_file( fn, mmap=False ).variables[ self.var ][ :, lev, ... ] for fn in self.filelist ])
	def get_dates( self ):
		dates = [ os.path.basename( fn )[:-3].split( '_' )[ len(os.path.basename( fn )[:-3].split( '_' )) - 1 ].split( '-' ) \
						for fn in self.filelist ]
		dates =[ j for i in dates for j in i ]
		dates = [ min( dates ), max( dates ) ]
		begin, end = [ dt.date( int(i[:4]), int(i[4:]), 15 ) for i in dates ]
		return [ i for i in rrule( MONTHLY, dtstart=begin, until=end ) ]
	def get_subset_data( self ):
		data = self.get_data( )
		dates = self.get_dates( )
		begin = dates[ 0 ]
		end = dates[ len( dates ) - 1 ]
		sub_begin = len( [ i for i in rrule( MONTHLY, dtstart=begin, until=self.begin_subset_date ) ] ) - 1
		sub_end = len( [ i for i in rrule( MONTHLY, dtstart=begin, until=self.end_subset_date ) ] )
		date_sub = dates[ sub_begin:sub_end ]
		data_sub = data[ sub_begin:sub_end, ... ]
		return date_sub, data_sub
	def calculate_climatology( self ):
		begin = len( [ i for i in rrule( MONTHLY, dtstart=self.begin_subset_date, until=self.climatology_begin_date ) ] ) - 1
		end = len( [ i for i in rrule( MONTHLY, dtstart=self.begin_subset_date, until=self.climatology_end_date ) ] )
		data_sub = self.data[ begin:end, ... ]
		clim_arr = np.dstack([ np.mean( data_sub[ j, ... ], axis=0 ) for j in [ range( i-1, 360, 12 ) for i in range( 1, 12+1, 1 ) ] ])
		clim_arr = np.rollaxis( clim_arr, axis=-1 )
		return clim_arr
	def calculate_anomalies( self ):
		return self.data - np.repeat( self.calculate_climatology(), self.data.shape[0]/12, axis=0 )



if __name__ == '__main__':

	import os, sys, re, glob, rasterio, fiona, shapely
	import numpy as np
	import pandas as pd
	import datetime as dt
	from dateutil.rrule import rrule, MONTHLY
	from scipy.io import netcdf_file
	import argparse

	# lets get some data to work with 
	# os.chdir( '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/cmip5/output1/IPSL/IPSL-CM5A-LR/historical/mon/atmos/Amon/r1i1p1/v20110406/hur' )
	# '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/cmip5/output1/IPSL/IPSL-CM5A-LR/historical/mon/atmos/Amon/r1i1p1/v20110406/hur',
	args = { 'input_dir': '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/cmip5/output1/NOAA-GFDL/GFDL-CM3/historical/mon/atmos/Amon/r1i1p1/v20120227/hur',
			 'level': 10000,
			 'var': 'hur',
			 'begin_subset_date':dt.date( 1900, 01, 15 ), 
			 'end_subset_date':dt.date( 2005, 12, 15 ),
			 'climatology_begin_date':dt.date( 1961, 01, 15 ),
			 'climatology_end_date':dt.date( 1990, 12, 15 ) }

	hur = RelativeHumidityDS( **args )


# here is some new stuff that might get us to the finish line with less code
from netCDF4 import Dataset, MFDataset, num2date
import pandas as pd
import numpy as np
import os, sys, re, xray


# input_dir = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/cmip5/output1/NOAA-GFDL/GFDL-CM3/historical/mon/atmos/Amon/r1i1p1/v20120227/hur'
input_dir = '/home/malindgren/Documents/hur'
# ds = MFDataset( os.path.join( input_dir, 'hur_*.nc' ) )
os.chdir( input_dir )

ds = Dataset( 'hur_Amon_GFDL-CM3_historical_r1i1p1_186001-186412.nc' )
atmos_level = 11

time = ds.variables[ 'time' ]
hur = ds.variables[ 'hur' ][ :, atmos_level, ... ]
lat = ds.variables[ 'lat' ][:]
lon = ds.variables[ 'lon' ][:]

# convert the dates in the file to python datetime
dates = num2date( time[:], time.units )
# convert the dates to a pandas datetime object from the base datetime objects generated from num2date
dates_pd = pd.to_datetime( dates )


# # # THIS IS THE NEW AND WORKING PORTION OF THIS SCRIPT

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

# # # # # # END KNOWN WORKING PORTION ## # # # # #


# we need to read in and reproject the cru data to match the desired output template
cru_data = '/workspace/Data/Base_Data/Climate/World/CRU_grids/CRU_TS20/grid_10min_reh.dat.gz'
output_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts20'

# read in the gzipped .dat file downloaded from the MET Office UK
cru_df = pd.read_csv( cru_data, delim_whitespace=True, compression='gzip', header=None, names=['lat', 'lon', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'] )
# lons = np.array( cru_df.lon.unique() )
# lons.sort()
# lats = cru_df.lat.unique()
# lats[::-1].sort()
# lats = np.array( lats )
# lats.sort()
# points = np.meshgrid( lons, lats )
# dat_cols = [u'01', u'02', u'03', u'04', u'05', u'06', u'07', u'08', u'09', u'10', u'11', u'12']
# from scipy.interpolate import RegularGridInterpolator

# xmin,xmax,ymin,ymax = [min(x),max(x),min(y),max(y)]

# #size of 1 m grid
# nx = (int(xmax - xmin + 1))#CHANGE HERE
# ny = (int(ymax - ymin + 1))#CHANGE HERE

# # Generate a regular grid to interpolate the data.
# xi = np.linspace(xmin, xmax, nx)
# yi = np.linspace(ymin, ymax, ny)
# xi, yi = np.meshgrid(xi, yi) 
 

# lonlat = cru_df[['lon', 'lat']]
# lonlat = lonlat.sort(['lon'])
# lons, lats = lonlat.iteritems()
# lons = np.array( lons[1] )
# lats = np.array( lats[1] )
# new = RegularGridInterpolator( ( lons, lats ), np.array( cru_df['01'] ), method='linear' )

# resolution = 0.16667
# xmin, xmax = lons.min(), lons.max()
# ymin, ymax = lats.min(), lats.max()

# this would be the entire world grid at 10' resolution
cols=2160
rows=1080
crs = 'epsg:4326'
resolution = 0.167
affine = A( *[resolution, 0.0, -180.0, 0.0, -resolution, 90.0] )
# count = time_len

meta = { 'affine':affine,
		'height':rows,
		'width':cols,
		'crs':crs,
		'driver':'GTiff',
		'dtype':np.float32,
		'count':1,
		'transform':affine.to_gdal(), 
		'compress':'lzw' }

# def world2Pixel(geoMatrix, x, y):
# 	"""
# 	Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
# 	the pixel location of a geospatial coordinate
# 	(source http://www2.geog.ucl.ac.uk/~plewis/geogg122/vectorMask.html)
# 	geoMatrix
# 	[0] = top left x (x Origin)
# 	[1] = w-e pixel resolution (pixel Width)
# 	[2] = rotation, 0 if image is "north up"
# 	[3] = top left y (y Origin)
# 	[4] = rotation, 0 if image is "north up"
# 	[5] = n-s pixel resolution (pixel Height)

# 	"""
# 	ulX = geoMatrix[0]
# 	ulY = geoMatrix[3]
# 	xDist = geoMatrix[1]
# 	yDist = geoMatrix[5]
# 	rtnX = geoMatrix[2]
# 	rtnY = geoMatrix[4]
# 	pixel = np.round((x - ulX) / xDist).astype(np.int)
# 	line = np.round((ulY - y) / xDist).astype(np.int)
# 	return (pixel, line)

# lats = cru_df.lat.tolist()
# lons = cru_df.lon.tolist()
# lonlat = zip(lons, lats)
# locations = [ world2Pixel( affine.to_gdal(), x, y ) for x, y in lonlat ]

# # break those locations into something that numpy can use
# xs = [x for x,y in locations ]
# ys = [y for x,y in locations ]
# locations = ( np.array( ys ), np.array( xs ) )

# # make a template and fill it!
# template = np.zeros( (rows, cols), dtype=np.float32 )
# template[ locations ] = np.array(cru_df['01']).tolist()

# # now just do this as a loop through the 12 months and output the files
# output_filename = '/home/UA/malindgren/Documents/test_out_cru_jan.tif'
# with rasterio.open( output_filename, 'w', **meta ) as out:
# 	out.write_band( 1, template.astype( np.float32 ) )

# # here we need to regrid to the AKCAN extent and resolution
# rst = rasterio.open( '/home/UA/malindgren/Documents/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif' )
# dst = np.empty_like( rst.read( 1 ) )

# reproject( template,
# 			dst,
# 			src_transform=affine.to_gdal(),
# 			src_crs={'init':'epsg:4326'},
# 			dst_transform=rst.affine.to_gdal(),
# 			dst_crs={'init':'epsg:3338'},
# 			resampling=RESAMPLING.cubic_spline )

# # # # # Lets try interp2d
# from scipy.interpolate import interp2d

# lats = cru_df.lat.tolist()
# lons = cru_df.lon.tolist()
# x, y = np.meshgrid(np.unique(lats), np.unique(lons))

# z = cru_data['01']

# interp2d( np.array( lats ), np.array( lons ), z, 'cubic' )




# # # # THIS IS THE WAY WE AER GOING TO DO THIS INTERPOLATION

from scipy.interpolate import griddata, Rbf
import matplotlib.mlab as mlab
import numpy as np

resolution = 0.16
x = np.array(cru_df.lon.tolist())
y = np.array(cru_df.lat.tolist())
z = np.array(cru_df['01'].tolist())
# lonlat = np.vstack( [x, y] ).T
lonlat = np.fliplr(np.array(cru_df[['lat','lon']]))

xi = np.linspace( x.min(), x.max(), 2160 )
yi = np.linspace( y.min(), y.max(), ny )
xi, yi = np.meshgrid( xi, yi )

# from scipy.interpolate import Rbf
# rbf = Rbf(x, y, z, epsilon=2)

# resolution = 0.16667
# xmin, xmax = min(lats) - resolution, max(lats) + resolution
# ymin, ymax = min(lons) - resolution, max(lons) + resolution

ymin,ymax,xmin,xmax = [ min(y), max(y), min(x), max(x) ]

resolution = 0.16666666666666666
x = np.array(cru_df.lon.tolist())
y = np.array(cru_df.lat.tolist())
z = np.array(cru_df['01'].tolist())
# lonlat = np.vstack( [x, y] ).T
lonlat = np.fliplr(np.array(cru_df[['lat','lon']]))

# grid size
# nx = 2159 #(int((xmax - xmin)/resolution) + 2) # CHANGE HERE
# ny = (int((ymax - ymin)/resolution) + 1) # CHANGE HERE

(ymax - ymin) / nrows
(xmax - xmin) / ncols

nrows = 900
ncols = 2160
ymin = -60.0
ymax = 90.0
xmin = -180.0
xmax = 180.0

# Generate a regular grid to interpolate the data.
xi = np.linspace( xmin, xmax, ncols )
yi = np.linspace( ymin, ymax, nrows )
xi, yi = np.meshgrid( xi, yi )


# cols=2160
# rows=860
# crs = 'epsg:4326'
# resolution = 0.167
affine = A( *[resolution, 0.0, xmin, 0.0, -resolution, ymax] )
# count = time_len

meta = { 'affine':affine,
		'height':nrows,
		'width':ncols,
		'crs':crs,
		'driver':'GTiff',
		'dtype':np.float32,
		'count':1,
		'transform':affine.to_gdal(), 
		'compress':'lzw' }

# def normalize_x(data):
#     data = data.astype(np.float)
#     return (data - xmin) / (xmax - xmin)

# def normalize_y(data):
#     data = data.astype(np.float)
#     return (data - ymin) / (ymax - ymin)

# x_new, xi_new = normalize_x(lats), normalize_x(xi)
# y_new, yi_new = normalize_y(lons), normalize_y(yi)

# lonlat_new = np.vstack( [x_new, y_new] ).T

# what we need to do here is:
# 1. reproject akcan extent to 4326
# 2. crop cru to that extent
# 3. project points to 3338
# 4. run the griddata on these outputs


# this is a snippet from an R script I found for working with this annoying dataset
wrld <- raster(nrows = 900, ncols = 2160, ymin = -60.0, ymax = 90.0, xmin = -180.0, xmax = 180.0)

zi = griddata( lonlat, np.array(cru_df['01'].tolist()), (xi,yi) , method='linear' )

output_filename = '/home/UA/malindgren/Documents/test_out_cru_jan_new2.tif'
with rasterio.open( output_filename, 'w', **meta ) as out:
	out.write_band( 1, np.flipud( zi.astype( np.float32 ) ) )


out = rasterio.open( out.name )

os.system( 'gdalwarp -of GTiff -co COMPRESS=LZW -wo SOURCE_EXTRA=1200 -overwrite -multi -tr 2000 2000 -te ' + '-2173223.206087799, 176412.93264414696, 4262776.793912201, 2548412.932644147' + ' -s_srs EPSG:4326 -t_srs EPSG:3338 ' + out.name + ' /home/UA/malindgren/Documents/test_gdalwarp2.tif' )

rst = rasterio.open( '/home/UA/malindgren/Documents/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif' )
dst = np.empty_like( rst.read( 1 ) )

reproject( out.read_band(1),
			dst,
			src_transform=out.affine.to_gdal(),
			src_crs={'init':'epsg:4326'},
			dst_transform=rst.affine.to_gdal(),
			dst_crs={'init':'epsg:3338'},
			resampling=RESAMPLING.bilinear,
			kwargs={'SOURCE_EXTRA':'120'} )

meta = rst.meta
meta.update( compress='lzw' )
output_filename = '/home/UA/malindgren/Documents/test_out_cru_jan_akcan.tif'
with rasterio.open( output_filename, 'w', **meta ) as out2:
	out2.write_band( 1, dst.astype( np.float32 ) )


# # # # RECTBIVARIATE SPLINE TEST
cru_df_sorted = cru_df.sort(['lat', 'lon'])

x = np.unique( np.array( cru_df_sorted.lon.tolist() ) )
y = np.unique( np.array( cru_df_sorted.lat.tolist() ) )
z = np.array(cru_df_sorted[ '01' ].tolist()).reshape( (nx, ny) )

ymin,ymax,xmin,xmax = [ min(y), max(y), min(x), max(x) ]

tmp = RectBivariateSpline(x, y, z, bbox=[-180.0, 180.0, -90.0, 90.0] )



# # # # # #
# Generate a regular grid to interpolate the data.
cols = 2160
rows = 1080
xi = np.linspace( -180.0, 180.0, cols )
yi = np.linspace( ymin, ymax, rows )
xi, yi = np.meshgrid( xi, yi )

# lats, lons = np.meshgrid( np.unique(lats), np.unique(lons) )
zi = SmoothBivariateSpline( x, y, np.array( cru_df['01'].tolist() ), kx=1, ky=1, bbox=[-180,180,-90,90] )
znew = zi.ev( xi.ravel(), yi.ravel() ).reshape( nx, ny )

output_filename = '/home/UA/malindgren/Documents/test_out_cru_jan_spline.tif'
with rasterio.open( output_filename, 'w', **meta ) as out:
	out.write_band( 1, np.flipud( znew.astype( np.float32 ) ) )


# Normalize coordinate system
def normalize_x(data):
    data = data.astype(np.float)
    return (data - xmin) / (xmax - xmin)

def normalize_y(data):
    data = data.astype(np.float)
    return (data - ymin) / (ymax - ymin)

x_new, xi_new = normalize_x(x), normalize_x(xi)
y_new, yi_new = normalize_y(y), normalize_y(yi)
# anisotropic
zi = mlab.griddata(x_new, y_new, z, xi_new, yi_new)
# this may do it isotropic
zi = mlab.griddata( lats, lons, np.array(z.tolist()), xi, yi,  )

plt.imshow( np.flipud( zi.T ) )
# plt.axis([xmin, xmax, ymin, ymax])
plt.colorbar()
plt.savefig( '/home/UA/malindgren/Documents/test_out13.png' )
plt.close()

zi = griddata( (lats, lons), np.array( cru_df_sorted['01'].tolist() ), (xi, yi), method='linear')

with rasterio.open( output_filename, 'w', **meta ) as out:
	out.write_band( 1, zi )

zi = SmoothBivariateSpline( lats, lons, np.array( cru_df['01'].tolist() ), bbox=[xmin,xmax,ymin,ymax] )
znew = zi.ev( xnew.ravel(), ynew.ravel() )
# this is another type of interp
import matplotlib.mlab as ml
zi = ml.griddata(lon,lat,jan,xi,yi,interp='nn') 

reproject( zi,
			dst,
			src_transform=affine.to_gdal(),
			src_crs={'init':'epsg:4326'},
			dst_transform=rst.affine.to_gdal(),
			dst_crs={'init':'epsg:3338'},
			resampling=RESAMPLING.cubic_spline )

output_filename = '/home/UA/malindgren/Documents/test_out_cru_jan3.tif'
with rasterio.open( output_filename, 'w', **meta ) as out:
	out.write_band( 1, np.flipud(zi.T).astype( np.float32 ) )

def cru_xyz_to_shp( in_xyz, lon_col, lat_col, crs, output_filename ):
	colnames = ['lat', 'lon', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
	from shapely.geometry import Point
	import pandas as pd
	import geopandas as gpd
	import os

	if os.path.splitext( in_xyz )[1] == '.gz':
		cru_df = pd.read_csv( cru_data, delim_whitespace=True, compression='gzip', header=None, names=colnames )
	else:
		cru_df = pd.read_csv( cru_data, delim_whitespace=True, header=None, names=colnames )

	# latlong = cru_df[ ['lon', 'lat'] ]
	# cru_df[ 'geometry' ] = [ Point( (j.lon, j.lat) ) for i,j in latlong.iterrows() ]

	# create a column named geometry with shapely geometry objects for each row
	def f( x ):
		''' return a Point shapely object for each x,y pair'''
		return Point( x.lon, x.lat )
	
	cru_df[ 'geometry' ] = cru_df.apply( f, axis=1 )

	# convert that to a GeoDataFrame object from GeoPANDAS
	cru_df = gpd.GeoDataFrame( cru_df )
	cru_df.to_file( output_filename, 'ESRI Shapefile' )
	return output_filename


output_filename = '/home/UA/malindgren/Documents/test_out_cru.shp'
cru_xyz_to_shp( cru_data, 'lon', 'lat', {'init':'epsg:4326'}, output_filename )

######################################################################################################################################################################
def cru_xyz_to_shp( in_xyz, lon_col, lat_col, crs, output_filename ):
	'''
	convert the cru cl2.0 1961-1990 Climatology data to a shapefile.

	This can handle the .dat format even if compressed with .gzip extension.

	in_xyz = path to the .dat or .dat.gz downloaded cru cl2.0 file from UK Met Office site
	lon_col = 
	lat_col = 
	crs = proj4string or epsg code
	output_filename = string path to the output filename to be created
	'''
	colnames = ['lat', 'lon', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
	from shapely.geometry import Point
	import pandas as pd
	import geopandas as gpd
	import os

	if os.path.splitext( in_xyz )[1] == '.gz':
		cru_df = pd.read_csv( cru_data, delim_whitespace=True, compression='gzip', header=None, names=colnames )
	else:
		cru_df = pd.read_csv( cru_data, delim_whitespace=True, header=None, names=colnames )

	# latlong = cru_df[ ['lon', 'lat'] ]
	# cru_df[ 'geometry' ] = [ Point( (j.lon, j.lat) ) for i,j in latlong.iterrows() ]

	# create a column named geometry with shapely geometry objects for each row
	def f( x ):
		''' return a Point shapely object for each x,y pair'''
		return Point( x.lon, x.lat )
	
	cru_df[ 'geometry' ] = cru_df.apply( f, axis=1 )
	
	# convert that to a GeoDataFrame object from GeoPANDAS
	cru_df = gpd.GeoDataFrame( cru_df )
	cru_df.to_file( output_filename, 'ESRI Shapefile' )
	return output_filename


output_filename = '/home/UA/malindgren/Documents/test_out_cru.shp'
cru_xyz_to_shp( cru_data, 'lon', 'lat', {'init':'epsg:4326'}, output_filename )

def bounds_to_extent( bounds ):
	'''
	take input rasterio bounds object and return an extent
	'''
	l,b,r,t = bounds
	return [ (l,b), (r,b), (r,t), (l,t), (l,b) ]
def extent_to_shapefile( extent, output_shapefile, proj4string ):
	''' convert an extent to a shapefile using its proj4string '''
	import geopandas as gpd
	from shapely.geometry import Polygon
	gpd.GeoDataFrame( {'extent_id':1, 'geometry':Polygon( extent )}, index=[1], crs=proj4string ).to_file( output_shapefile, 'ESRI Shapefile' )
	return output_shapefile
def pad_bounds( rst, npixels, crs, output_shapefile ):
	'''
	convert the extents of 2 overlapping rasters to a shapefile with 
	an expansion of the intersection of the rasters extents by npixels

	rst1: rasterio raster object
	rst2: rasterio raster object
	npixels: tuple of 4 (left(-),bottom(-),right(+),top(+)) number of pixels to 
		expand in each direction. for 5 pixels in each direction it would look like 
		this: (-5. -5. 5, 5) or just in the right and top directions like this:
		(0,0,5,5).
	crs: epsg code or proj4string defining the geospatial reference 
		system
	output_shapefile: string full path to the newly created output shapefile

	'''
	import rasterio, os, sys
	from shapely.geometry import Polygon

	resolution = rst.res[0]
	new_bounds = [ bound+(expand*resolution) for bound, expand in zip( rst.bounds, npixels ) ]
	
	new_ext = bounds_to_extent( new_bounds )
	return extent_to_shapefile( new_ext, output_shapefile, crs )

# try to pad the bounds of the akcan extent to use in cropping and clipping the points using ogr2ogr
crs = { 'init':'epsg:3338' }
output_shapefile = '/home/UA/malindgren/Documents/akcan_extent.shp'
npixels = ( -200, -2000, 200, 200 )
pad_bounds( rst, npixels, crs, output_shapefile )

# # here we are going to reproject and crop to the AKCAN extent, the crudata built above
# extent_to_shapefile( bounds_to_extent( rst.bounds ), '/home/UA/malindgren/Documents/akcan_extent.shp', 'EPSG:3338') # this is from NSIDC HSIA stuff
# os.system( "ogr2ogr -wrapdateline -overwrite -f 'ESRI Shapefile' -clipdst '/home/UA/malindgren/Documents/akcan_extent.shp' -t_srs 'EPSG:4326' /home/UA/malindgren/Documents/test_out_cru_akcan.shp /home/UA/malindgren/Documents/test_out_cru.shp" )
os.system( "ogr2ogr -wrapdateline -overwrite -f 'ESRI Shapefile' -clipdst /home/UA/malindgren/Documents/akcan_extent.shp -s_srs 'EPSG:4326' -t_srs 'EPSG:3338' /home/UA/malindgren/Documents/test_out_cru_3338_expand.shp /home/UA/malindgren/Documents/test_out_cru.shp" )

rst = rasterio.open( '/home/UA/malindgren/Documents/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif' )
resolution = rst.res
meta = rst.meta

padded_bounds = fiona.open( output_shapefile ).bounds
xmin, ymin, xmax, ymax = padded_bounds
cols = (xmax - xmin) / resolution[1]
rows = (ymax - ymin) / resolution[0]

# update vars
meta[ 'affine' ] = A( resolution[0], 0.0, xmin, 0.0, -resolution[1], ymax )
meta[ 'crs' ] = { 'init':'epsg:3338' }
meta[ 'height' ] = rows
meta[ 'width' ] = cols
meta[ 'transform' ] = meta[ 'affine' ].to_gdal()

# dst = np.ndarray( (rows, cols) )

# read in the clipped and warped file from above
cru_df_akcan = gpd.read_file( '/home/UA/malindgren/Documents/test_out_cru_3338_expand.shp' )

# update lon and lat to the 3338
cru_df_akcan.lon = cru_df_akcan.geometry.apply( lambda x: x.x )
cru_df_akcan.lat = cru_df_akcan.geometry.apply( lambda x: x.y )

# build the interpolation input values
x = np.array(cru_df_akcan.lon.tolist())
y = np.array(cru_df_akcan.lat.tolist())
z = np.array(cru_df_akcan['01'].tolist())
# lonlat = np.vstack( [x, y] ).T
# lonlat = np.fliplr(np.array(cru_df[['lat','lon']]))

# xmin, ymin, xmax, ymax = rst.bounds
# rows, cols = rst.shape
xi = np.linspace( xmin, xmax, cols )
yi = np.linspace( ymin, ymax, rows )
xi, yi = np.meshgrid( xi, yi )

# run interpolation
zi = griddata( (x,y), np.array( cru_df_akcan['01'].tolist() ), (xi, yi), method='cubic' )
zi = np.flipud( zi.astype( np.float32 ) )

meta.update( compress='lzw' )
output_filename = '/home/UA/malindgren/Documents/hur_cru_cl20_jan.tif'
with rasterio.open( output_filename, 'w', **meta ) as out:
	out.write_band( 1, zi )

# what we really need to do here is not reproject, but crop down the data to the extent we need
# gdalwarp may be most appropriate? though there are ways to do with rasterio.
# see this post: https://github.com/mapbox/rasterio/issues/78

def window(self, *bbox):
	''' window from bounds -- sean gillies '''
    llidx = self.index(*bbox[:2])
    uridx = self.index(*bbox[2:])
    return (uridx[0], llidx[0]), (llidx[1], uridx[1])

# what needs to be developed here:
# 1. rasterio akcan bounds to window of the larger raster
window = bounds_to_window( rst.bounds )
out_arr = out.read_band( 1, window=window )

arr = rst.read( 1 )
mask = np.zeros_like( arr )
mask[ arr >= -90 ] = 1

# 2. mask the raster with the AKCAN mask
out_arr[ mask != 1 ] = -9999

# 3. build out a new meta for it based on the akcan raster
meta = rst.meta
meta.update( compress='lzw', crs=crs, nodata=-9999 )

# 4. write out the new cropped / masked raster
output_filename = '/home/UA/malindgren/Documents/hur_cru_cl20_jan_akcan.tif'
with rasterio.open( output_filename, 'w', **meta ) as out2:
	out2.write_band( 1, out_arr.astype( np.float32 ) )



# dst = np.empty_like( rst.read( 1 ) )
# reproject( zi,
# 			dst, 
# 			src_transform=out.affine.to_gdal(),
# 			src_crs={'init':'epsg:3338'},
# 			dst_transform=rst.affine.to_gdal(),
# 			dst_crs={'init':'epsg:3338'},
# 			resampling=RESAMPLING.cubic_spline )

# # write it out to a tif
# meta = rst.meta
# meta.update( compress='lzw', crs=crs )
# output_filename = '/home/UA/malindgren/Documents/hur_cru_cl20_jan_akcan.tif'
# with rasterio.open( output_filename, 'w', **meta ) as out2:
# 	out2.write_band( 1, dst.astype( np.float32 ) )


# template output dataset

# now is it possible to read this back in and interpolate on the na's?
# or should i be thinking of masking here and finding locations where the mask is incorrect 
#  and filling only those with interpolated values?  This would fill what we need and solve 
#  these current data woes...
# zi = griddata( (x,y), np.array( cru_df_akcan['01'].tolist() ), (xi, yi), method='linear' )



######################################################################################################################################################################




def hur2vap( hur_arr, tas_arr ):
	'''
	take a relative humidity 2d array and a temperature at 
	surface 2d array (same time slice) and compute vapor pressure

	returns:
		2d array of vapor pressure values across the domain

	NOTE: arrays must be the same shape and be datatypes which can
		properly operate on each other.

	'''
	esa = 6.112 * exp( (17.62*tas_arr) / (243.12+tas_arr) )
	return ( hur.mat * esa ) / 100




# loop through with a list comprehension the ndarray and apply the:
# reprojection
# ? units conversion ?
# downscaling
# output raster to disk (GTiff)

# return both the relative humidity and the vapor pressure for review? meh.
#  just return the vapor pressure. 
def hur_downscale( hur_arr, tas_arr, cru_arr, input_affine, output_affine, input_crs={'init':'epsg:4326'}, output_crs={'init':'epsg:3338'}, output_filename ):
	'''
	hur_arr and tas_arr are both 2d and must be the same shape and comparable datatypes
	they should both be anomalies as well (calculated prior).

	cru_arr is the corresponding higher-res climatology (in this case 10' TS2.0) which
	has been previously reprojected to the desired extent / resolution / shape.
	'''
	# this could be an argument
	resampler = RESAMPLING.cubic_spline
	# reproject to AKCAN 2km EPSG:3338
	# hur
	hur_akcan = np.empty_like( cru_arr )
	reproject( hur_arr,
				hur_akcan,
				src_transform=input_affine.to_gdal(),
				src_crs=input_crs,
				dst_transform=output_affine.to_gdal(),
				dst_crs=output_crs,
				resampling=resampler )
	# tas
	tas_akcan = np.empty_like( cru_arr )
	reproject( tas_arr,
				tas_akcan,
				src_transform=input_affine.to_gdal(),
				src_crs={'init':'epsg:4326'},
				dst_transform=output_affine.to_gdal(),
				dst_crs={'init':'epsg:3338'},
				resampling=resampler )

	del hur_arr, tas_arr # remove the input arrays

	# units conversion -> saturated vapor pressure
	vap_akcan = hur2vap( hur_akcan, tas_akcan )
	# now actually downscale it
	vap_akcan = cru_arr + vap_akcan

	with rasterio.open( output_filename, 'w', **meta ) as out_arr:
		out_arr.write( 1, vap_akcan )
	return output_filename







def calc_resolution( ds ):
	time_len, rows, cols = ds.shape


	latmin = min( ds.lat )
	latmax = max(ds.lat )
	lonmin = min( ds.lon )
	lonmax = max( ds.lon )

	(latmax - latmin) / rows
	(lonmax - lonmin)


# convert the sliced data to a pandas dataframe
df = hur_lev.to_dataframe()

# this is a simple way to rotate the coordinates of longitudes from 0-360 to -180-180
new_lons = ((ds.lon.data + 180) % 360) - 180

lons_greenwich = ((lons_pacific + 180) % 360) - 180

# see also:
# basemap shiftdata
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
m = Basemap( epsg=4326 )
new_arr = np.rollaxis( hur_lev.values.T, -1 )
test = m.shiftdata( xds.lon.data, hur_lev.values[1,...][:].T, lon_0=0.0 )

llons, llats = np.meshgrid( test[0], ds.lat.data )
x, y = m(llons, llats)
# m.contourf(x, y, test[1])
# m.drawcoastlines()
m.imshow( test[1] )
plt.savefig('test.png')


scp malindgren@137.229.94.86:/home/malindgren/Documents/hur/test.png ~/Documents

# from frank
gdalwarp -multi -overwrite -of NetCDF -t_srs WGS84 hur_Amon_GFDL-CM3_historical_r1i1p1_186001-186412.nc test.nc -wo SOURCE_EXTRA=1000 -wo NUM_THREADS=ALL_CPUS --config CENTER_LONG 0

	# test file
	# nc = netcdf_file( 'hur_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001-194912.nc' )
	# # get the pressure levels so we can query it for the one we want
	# # plev = nc.variables[ 'plev' ].data
	# lev, = np.where( nc.variables[ 'plev' ].data == 10000 )
	# lev = lev.tolist()[ 0 ] # grab that elements location as a single integer

	# # slice the 4d array to the 3dimensions we actually want for this exercise.
	# hur = nc.variables[ 'hur' ][ :, lev, ... ]

	# # now get the begin and end dates from the filename as dictionaries
	# time = nc.filename[ :-3 ].split( '_' )
	# time = time[ len( time ) - 1 ]
	# time = time.split( '-' )
	# begin, end = [ dt.date( int(i[:4]), int(i[4:]), 15 ) for i in time ]
	# all_dates = [ i for i in rrule( MONTHLY, dtstart=begin, until=end ) ]

	# # this is how we can get some new date information in string form
	# new_datefmt = i.strftime( "%m_%Y" )

rst = rasterio.open('/Data/Base_Data/Climate/AK_CAN_2km/projected/AR5_CMIP5_models/rcp26/CCSM4/tas/tas_mean_C_AR5_CCSM4_rcp26_01_2039.tif')
ext = '%d %d %d %d' % rst.bounds

fn1 = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/cmip5/output1/IPSL/IPSL-CM5A-LR/historical/mon/atmos/Amon/r1i1p1/v20110406/hur/hur_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001-194912.nc'
output_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/working'
os.system( 'gdalwarp -of GTiff -multi -overwrite -co NUM_THREADS=ALL -co COMPRESS=LZW -r bilinear -t_srs EPSG:3338 -tr 2000 2000 -tap -2173223 176412 4262776 2548412 ' + fn1 + ' ' + os.path.join( output_path, os.path.basename(fn1)))

# we need to use the python library to do the slicing of the needed pressure level for us to pass to gdalwarp and see if we can overcome a massive amount of time-consuming python
os.system( 'gdalwarp -of netCDF -r bilinear -t_srs EPSG:3338 -tr 2000 2000 -tap -te -2173223 176412 4262776 2548412 ' + 'NETCDF:'+fn1+':hur' + ' ' + os.path.join( output_path, os.path.basename(fn1) ) )


