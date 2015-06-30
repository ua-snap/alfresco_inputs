def read_ar5_mon( path, fn_prefix_filter, variable, level=None, time_begin='1900-01-01', time_end='2005-12-31', **kwargs ):
	'''
	open file(s) of a given netcdf dataset
	'''
	xds = xray.open_mfdataset( os.path.join( path, fn_prefix_filter ) )
	var = xds[ variable ].loc[ time_begin:time_end ]
	if level:
		out = var[ :, level, ... ]
	else:
		out = var[ :, ... ]
	return out

if __name__ == '__main__':
	import pandas as pd
	import numpy as np
	import os, sys, re, xray, rasterio
	from rasterio import Affine as A
	from rasterio.warp import reproject, RESAMPLING
	from mpl_toolkits.basemap import shiftgrid, addcyclic
	from functools import partial

	# NOTE: xray is currently only working with the pip install git+https://github.com/xray/xray

	# some setup pathing
	input_dir = '/home/UA/malindgren/Documents/hur'
	os.chdir( input_dir )
	path = input_dir
	fn_prefix_filter ='hur_Amon_GFDL-CM3_historical_r1i1p1_*.nc'
	variable = 'hur'
	time_begin='1900-01-01'
	time_end='2005-12-31'
	template_fn = '/home/UA/malindgren/Documents/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
	atmos_level = 11

	# open the data and subset to the needed atmos level
	ds = read_ar5_mon( path, fn_prefix_filter, variable, level=atmos_level )
	
	# generate climatology / anomalies
	climatology = ds.loc[ '1961-01-01':'1990-12-31' ].groupby( 'time.month' ).mean( 'time' )
	anomalies = ds.groupby( 'time.month' ) - climatology

	# some setup of the output raster metadata
	time_len, rows, cols = anomalies.shape
	crs = 'epsg:4326'
	affine = A( *[np.diff( ds.lon )[ 0 ], 0.0, -180.0, 0.0, -np.diff( ds.lat )[ 0 ], 90.0] )
	count = time_len
	resolution = ( np.diff( ds.lat )[ 0 ], np.diff( ds.lon )[ 0 ] )

	# shift the grid to Greenwich Centering
	dat, lons = shiftgrid( 180., anomalies[:], anomalies.lon.data, start=False )#

	# metadata for input?
	meta_4326 = { 'affine':affine,
			'height':rows,
			'width':cols,
			'crs':crs,
			'driver':'GTiff',
			'dtype':np.float32,
			'count':time_len,
			'transform':affine.to_gdal(),
			'compress':'lzw' }

	# grab template metadata and write and reproject to it
	template = rasterio.open( template_fn )
	meta_3338 = template.meta
	meta_3338.update( compress='lzw', crs={'init':'epsg:3338'} ) # 	meta_3338 = rasterio.open( template_fn ).meta

	# output_filename = 'test.tif'

	# with rasterio.open( output_filename, 'w', **meta_3338 ) as out:
	# 	print 'generated new output rst %s' % output_filename

	# out = []
	# for x in np.vsplit( dat, time_len ):
	# 	dst = np.empty_like( template.read(1) ) # np.zeros( shape=( time_len, meta_3338[ 'height' ], meta_3338[ 'width' ]))
	# 	reproject( x,
	# 				dst,
	# 				src_transform=meta_4326[ 'affine' ],
	# 				src_crs=meta_4326[ 'crs' ],
	# 				dst_transform=meta_3338[ 'affine' ],
	# 				dst_crs={'init':'epsg:3338'},
	# 				resampling=RESAMPLING.cubic_spline )
	# 	out.append( dst )

def downscale( src, dst, cru, src_crs, src_affine, dst_crs, dst_affine, output_filename, dst_meta, \
		method=RESAMPLING.cubic_spline, operation='add' **kwargs ):
	'''
	operation can be one of two keywords for the operation to perform the delta downscaling
	- keyword strings are one of: 'add'= addition, 'mult'=multiplication, or 'div'=division
	'''
	out = dst[:]
	# reproject src to dst
	reproject( src,
			out,
			src_transform=src_affine,
			src_crs=src_crs,
			dst_transform=dst_affine,
			dst_crs=dst_crs,
			resampling=method )
	
	{ 'add':np.sum, 'mult':np.mult, }

	if operation == 'add':
		downscaled = out + 



def f( src_hur, src_tas, src_cru, dst, src_crs, src_affine, dst_crs, dst_affine, output_filename, dst_meta, \
		method=RESAMPLING.cubic_spline, **kwargs ):
	
	# reproject / resample to akcan
	# hur
	hur = dst[:]
	reproject( src_hur,
				hur,
				src_transform=src_affine,
				src_crs=src_crs,
				dst_transform=dst_affine,
				dst_crs=dst_crs,
				resampling=method )

	# tas
	hur = dst[:]
	reproject( src_hur,
				hur,
				src_transform=src_affine,
				src_crs=src_crs,
				dst_transform=dst_affine,
				dst_crs=dst_crs,
				resampling=method )

	# downscale - add climatology and anomaly files for the same month
	# for temperature
	downscaled.month.tas <- cru.tas.current+nc.interp.r
	
	#### END TAS

	# # here we get the values of the current downscaled layer
	downscaled.month.flip.hur.r.v <- getValues(downscaled.month.flip.hur.r)
	
	# now we change the >100 values to 95 as per steph's okaying it
	values(downscaled.month.flip.hur.r)[which(values(downscaled.month.flip.hur.r) > 100)] <- 95

	print("***** 	COVNERTING RELATIVE HUMIDITY BACK TO VAPOR PRESSURE 	******************")

	# make the rasters into matrices
	tas.mat <- getValues(downscaled.month.flip.tas.r) 
	hur.mat <- getValues(downscaled.month.flip.hur.r)

	# new rasters to hold the output based on the information from the current set
	vapor.ts31.10min <- raster(downscaled.month.flip.tas.r)

	# convert back to vapor pressure
	esa = 6.112*exp(17.62*tas.mat/(243.12+tas.mat))
	vapor.ts31.10min.hold = (hur.mat*esa)/100
	# set the converted values into a vapor pressure variable
	vapor.ts31.10min <- setValues(vapor.ts31.10min,vapor.ts31.10min.hold)

	with rasterio.open( output_filename, 'w', **dst_meta ) as out:
		out.write_band( 1, dst2 )


	return output_filename

dst = np.empty_like( rasterio.open( template_fn ).read( 1 ) )
# from pathos import multiprocessing as mp
import multiprocessing as mp
pool = mp.Pool( 10 )
# f2 = partial( f, dst=dst, src_crs=meta_4326[ 'crs' ], src_affine=meta_4326[ 'affine' ], dst_crs=meta_3338[ 'crs' ], dst_affine=meta_3338[ 'affine' ], dst_meta=meta_3338 )
output_dir = '/home/UA/malindgren/Documents/hur/akcan'
output_filenames = [ os.path.join( output_dir, 'hur_level11_akcan_'+str(i) + '.tif' ) for i in range( time_len ) ]
out = pool.map( lambda x: f( **x ), [{'src':src, 'output_filename':fn, 'dst':dst, 'src_crs':meta_4326[ 'crs' ], 'src_affine':meta_4326[ 'affine' ], 'dst_crs':meta_3338[ 'crs' ], 'dst_affine':meta_3338[ 'affine' ], 'dst_meta':meta_3338 } for src,fn in zip( np.vsplit( dat, time_len ), output_filenames ) ] )
pool.close()

# # # 

# open multiple datasets as a single file
xds = xray.open_mfdataset( 'hur_Amon_GFDL-CM3_historical_r1i1p1_*.nc' )
xds_hur = xds.hur.loc[ '1900-01-01':'2005-12-12' ] # slice the dataset using the time variable in xray object
hur_lev = xds_hur[ :, atmos_level, ... ]

# calculate climatology and anomalies
climatology = hur_lev.loc[ '1961-01-01':'1990-12-31' ].groupby( 'time.month' ).mean( 'time' )
anomalies = hur_lev.groupby( 'time.month' ) - climatology

# # # REPROJECT AND CROP EXTENT

time_len, rows, cols = hur_lev.shape
# NOTE: geotransform = [left, res, 0.0, top, 0.0, res]
height = rows
width = cols
crs = 'epsg:4326'
affine = A( *[np.diff( xds.lon )[ 0 ], 0.0, -180.0, 0.0, -np.diff( xds.lat )[ 0 ], 90.0] )
count = time_len
resolution = ( np.diff( xds.lat )[ 0 ], np.diff( xds.lon )[ 0 ] )

# shift the grid to Greenwich Centering
dat, lons = shiftgrid( 180., anomalies[:], anomalies.lon.data, start=False )#

# metadata for input?
meta = { 'affine':affine,
		'height':rows,
		'width':cols,
		'crs':crs,
		'driver':'GTiff',
		'dtype':np.float32,
		'count':1 }

meta.update( transform=meta[ 'affine' ].to_gdal() )

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

# facets
base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data'
project = 'cmip5'
institute = 'IPSL' # i wonder if this one is really needed or can we just traverse past it?
model = 'IPSL-CM5A-LR'  # this is an upper-case representation and can
experiment = 'historical'
frequency = 'mon'
realm = 'atmos'
cmor_table = 'Amon'
ensemble = 'r1i1p1'
variable = 'hur'

hur_Amon_rcp45_r1i1p1_v20110914

patterns = ['*'.join([project, institute, model, experiment, frequency, realm, cmor_table, ensemble, variable])]
find_files( base_path, patterns )

# /cmip5/output1/IPSL/IPSL-CM5A-LR/historical/mon/atmos/Amon/r1i1p1/v20110406

def rec_dd():
    return defaultdict(rec_dd)

import os
from collections import defaultdict

dd = defaultdict( rec_dd )
for root, subdir, files in os.walk( base_path ):
	if len(files) is not 0:
		{'root':root, 'subdir':subdir, 'files':files }
		print [ os.path.join( root, subdir, f ) for f in files ]

# directory traversal
import fnmatch
import functools
import itertools
import os

# Remove the annotations if you're not on Python3
def find_files( dir_path, patterns ):
    """
    Returns a generator yielding files matching the given patterns
    :type dir_path: str
    :type patterns: [str]
    :rtype : [str]
    :param dir_path: Directory to search for files/directories under. Defaults to current dir.
    :param patterns: Patterns of files to search for. Defaults to ["*"]. Example: ["*.json", "*.xml"]
    """
    # not sure why these exist...
    path = dir_path
    if not patterns:
    	path_patterns = [ "*" ]
    else:
	    path_patterns = patterns

    for root_dir, dir_names, file_names in os.walk(path):
        filter_partial = functools.partial(fnmatch.filter, file_names)

        for file_name in itertools.chain( *map( filter_partial, path_patterns ) ):
            yield os.path.join( root_dir, file_name )

# groupby something in the new series
def grouping_files( x ):
	dir_path = os.path.dirname( x )
	fn, _ = os.path.splitext( os.path.basename( x ) )
	# remove dates from filename -- they have a hyphen
	fn_base = '_'.join([ i for i in fn.split( '_' ) if '-' not in i ])
	# return the path element that startswith 'v' this is the version attribute
	version = [ x for x in dir_path.split( os.path.sep ) if x.startswith( 'v' ) ]
	return '_'.join([ fn_base, version[0] ])

def group_versions( x ):
	dir_path = os.path.dirname( x )
	# return the path element that startswith 'v' this is the version attribute
	version = [ x for x in dir_path.split( os.path.sep ) if x.startswith( 'v' ) ]
	return version

matches = pd.Series([ match for match in find_files( base_path, [ '*hur_*' )] ])

hold = [ group for group in matches.groupby( matches.apply( grouping_files ) )]

hold = [ group for group in matches.groupby( matches.apply( grouping_files ) ).groupby( matches.apply( group_versions ) )]

# OR do we just do the simple task of only returning a single 
# https://github.com/tsileo/dirtools look at this package


# 2 stage groupings:
# get all the data for a particular variable
# group all data by the model name
# group each group by the version number





