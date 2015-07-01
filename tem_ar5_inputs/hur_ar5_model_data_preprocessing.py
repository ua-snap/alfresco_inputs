def read_ar5_mon( path, fn_prefix_filter, variable, level=None, time_begin='1900-01-01', time_end='2005-12-31', **kwargs ):
	'''
	open file(s) of a given netcdf dataset
	using xray mfdataset
	'''
	xds = xray.open_mfdataset( os.path.join( path, fn_prefix_filter ) )
	var = xds[ variable ].loc[ time_begin:time_end ]
	if level:
		out = var[ :, level, ... ]
	else:
		out = var[ :, ... ]
	return out
def cru_generator( n, cru_clim_list ):
	'''
	generator that will produce the cru climatologies with a
	generator and replicate for the total number of years in n
	'''
	months = [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12' ]
	for i in range( n ):
		for count, j in enumerate( cru_clim_list ):
			yield j
def downscale( src, dst, cru, src_crs, src_affine, dst_crs, dst_affine, output_filename, dst_meta, \
		method=RESAMPLING.cubic_spline, operation='add', **kwargs ):
	'''
	operation can be one of two keywords for the operation to perform the delta downscaling
	- keyword strings are one of: 'add'= addition, 'mult'=multiplication, or 'div'=division (not implemented)
	'''
	def add( cru, anom ):
		return cru + anom		
	def mult( cru, anom ):
		return cru * anom
	def div( cru, anom ):
		# return cru / anom
		# this one may not be useful, but the placeholder is here 
		return NotImplementedError

	# reproject src to dst
	out = dst[:]
	reproject( src,
			out,
			src_transform=src_affine,
			src_crs=src_crs,
			dst_transform=dst_affine,
			dst_crs=dst_crs,
			resampling=method )
	# switch to deal with different downscaling operators
	operation_switch = { 'add':add, 'mult':mult, 'div':div }
	downscaled = operation_switch[ operation ]( cru, out )

	# this is a geotiff creator so lets pass in the lzw compression
	dst_meta.update( compress='lzw' )
	with rasterio.open( output_filename, 'w', **dst_meta ) as out:
		out.write_band( 1, downscaled )
	return output_filename
def run( args ):
	''' 
	simple function wrapper for unpacking an argument dict 
	to the downscale function for getting around the single 
	argument pass to multiprocessing.map implementation issue.
	'''
	return( downscale( **args ) )

if __name__ == '__main__':
	import pandas as pd
	import numpy as np
	import os, sys, re, xray, rasterio, glob
	from rasterio import Affine as A
	from rasterio.warp import reproject, RESAMPLING
	from mpl_toolkits.basemap import shiftgrid, addcyclic
	from pathos import multiprocessing as mp

	# NOTE: xray is currently only working with: pip install git+https://github.com/xray/xray
	#		install pathos with: pip install git+https://github.com/uqfoundation/pathos

	# some setup pathing
	input_dir = '/home/UA/malindgren/Documents/hur'
	os.chdir( input_dir )
	path = input_dir
	fn_prefix_filter ='hur_Amon_GFDL-CM3_historical_r1i1p1_*.nc'
	variable = 'hur'
	time_begin='1900-01-01'
	time_end='2005-12-31'
	template_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
	atmos_level = 11
	cru_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts20/akcan'

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
	meta_4326 = {'affine':affine,
				'height':rows,
				'width':cols,
				'crs':crs,
				'driver':'GTiff',
				'dtype':np.float32,
				'count':time_len,
				'transform':affine.to_gdal(),
				'compress':'lzw' }

	# grab template metadata and write and reproject to it
	meta_3338 = rasterio.open( template_fn ).meta
	meta_3338.update( compress='lzw', crs={'init':'epsg:3338'} )

	dst = np.empty_like( rasterio.open( template_fn ).read( 1 ) )
	output_dir = '/home/UA/malindgren/Documents/hur/akcan/new'
	output_filenames = [ os.path.join( output_dir, 'hur_level11_akcan_' + str(i) + '.tif' ) for i in range( time_len ) ]
	cru_files = glob.glob( os.path.join( cru_path, '*.tif' ) )
	cru_files.sort()

	# this is effectively working
	cru_stack = [ rasterio.open( fn ).read( 1 ) for fn in cru_files ]
	cru_gen = cru_generator( 1271, cru_stack )

	# run in parallel
	pool = mp.Pool( 10 )
	out = pool.map( run, [{'src':src, 'output_filename':fn, 'dst':dst, 'cru':cru, 'src_crs':meta_4326[ 'crs' ], 'src_affine':meta_4326[ 'affine' ], 'dst_crs':meta_3338[ 'crs' ], 'dst_affine':meta_3338[ 'affine' ], 'dst_meta':meta_3338, 'operation':'add' } for src,fn,cru in zip( np.vsplit( dat, time_len ), output_filenames, cru_gen ) ] ) 
	pool.close()



# this should be integrated at the top-level
# directory traversal <<- this stuff is now working 
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
	import itertools, functools
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

# lets try with models
models = [ 'GISS-E2-R',  ] # fill this
model = 'GISS-E2-R'

# get all matches first
matches = pd.Series([ match for match in find_files( base_path, [ 'hur_*'+model+'*' ] ) ])

# now lets group em by version
grouped = dict([ group for group in matches.groupby( matches.apply( grouping_files ) )])

# now we need the keys so that we can split the strings into attributes to parse
keys_df = pd.DataFrame({ key:key.split( '_' ) for key in grouped.keys() }).T

def drop_old_versions( df ):
	rows,cols = df.shape
	if rows > 1 & rows < 3:
		version_nums = df[ 4 ].apply( lambda x : int( x.replace( 'v', '' ) ) )
		# max( version_nums )
		return df.drop( df[df[4] != 'v' + str( max( version_nums ) )].index )
	elif rows > 3:
		# potentially unnecessary
		None
	else:
		return df

# parse the keys values to keep only the ones that are latest versions
keys_df_grouped = pd.concat([ drop_old_versions(i[1]) for i in keys_df.groupby( 2 ) ])

# now keep only the keys we want
final_out = { k:v for k,v in grouped.iteritems() if k in keys_df_grouped.index.tolist() }



# # # # THIS IS A TESTING AREA TO FIGURE OUT THE BEST WAY TO PRESENT THE ALGORITHM WITH DATA FROM THE HOLDINGS
# THIS IS NOT YET COMPLETE!

# facets
# base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data'
# project = 'cmip5'
# institute = 'IPSL' # i wonder if this one is really needed or can we just traverse past it?
# model = 'IPSL-CM5A-LR'  # this is an upper-case representation and can
# experiment = 'historical'
# frequency = 'mon'
# realm = 'atmos'
# cmor_table = 'Amon'
# ensemble = 'r1i1p1'
# variable = 'hur'

# patterns = ['*'.join([project, institute, model, experiment, frequency, realm, cmor_table, ensemble, variable])]
# find_files( base_path, patterns )

# /cmip5/output1/IPSL/IPSL-CM5A-LR/historical/mon/atmos/Amon/r1i1p1/v20110406

# def rec_dd():
# 	return defaultdict(rec_dd)

# import os
# from collections import defaultdict

# dd = defaultdict( rec_dd )
# for root, subdir, files in os.walk( base_path ):
# 	if len(files) is not 0:
# 		{'root':root, 'subdir':subdir, 'files':files }
# 		print [ os.path.join( root, subdir, f ) for f in files ]


