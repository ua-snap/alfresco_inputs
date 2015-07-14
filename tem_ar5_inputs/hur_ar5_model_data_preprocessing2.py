# # # # #
# Tool to downscale the CMIP5 data from the PCMDI group. 
# # # # #

def read_ar5_mon( path, variable, filelist=None, fn_prefix_filter=None, level=None, time_begin='1900-01-01', time_end='2005-12-31', **kwargs ):
	'''
	open file(s) of a given netcdf dataset
	using xray mfdataset
	'''
	if filelist and fn_prefix_filter:
		NameError( 'both filelist and fn_prefix_filter are included, only one can be used at time' )
	elif filelist:
		xds = xray.open_mfdataset( filelist )
	elif fn_prefix_filter:
		xds = xray.open_mfdataset( os.path.join( path, fn_prefix_filter ) )
	else:
		NameError( 'one of filelist or fn_prefix_filter must be included' )
	
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
	import os, sys, re, xray, rasterio, glob, argparse
	from rasterio import Affine as A
	from rasterio.warp import reproject, RESAMPLING
	from mpl_toolkits.basemap import shiftgrid, addcyclic
	from pathos import multiprocessing as mp

	# NOTE: xray is currently only working with: pip install git+https://github.com/xray/xray
	#		install pathos with: pip install git+https://github.com/uqfoundation/pathos

	# some setup pathing <<-- THIS TO BE CONVERTED TO ARGUMENTS AT COMMAND LINE
	input_dir = '/home/UA/malindgren/Documents/hur'
	output_dir = '/home/UA/malindgren/Documents/hur/akcan/new'
	path = input_dir
	fn_prefix_filter ='hur_Amon_GFDL-CM3_historical_r1i1p1_*.nc'
	variable = 'hur'
	time_begin = '1900-01-01' # will change for future and historical
	time_end = '2005-12-31' # will change for future and historical
	climatology_begin = '1961-01-01'
	climatology_end = '1990-12-31'
	template_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
	atmos_level = 11
	cru_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts20/akcan'

	# upack argparse'd inputs from the command line


	# A VERY IMPORTANT NOTE HERE IS THAT WE ARE GOING TO ALSO HAVE TO READ IN THE HISTORICAL DATA FOR 
	# PERIOD 1961-1990 BECAUSE WE NEED TO CALCULATE ANOMALIES FROM THAT WITH THE MODELED FUTURES.
	# THIS IS GOING TO REQUIRE SOME THOUGHT ABOUT HOW WE WILL INTEGRATE THIS INTO A SINGLE, DYNAMIC SCRIPT.

	os.chdir( input_dir ) # maybe this should just come from the file list? dirname?

	# open the data and subset to the needed atmos level
	
	# THIS NEEDS TO BE MADE ABLE TO DEAL WITH AN INPUT OF TYPE LIST
	# WHERE LIST IS A GROUP OF FILENAMES TO BE USED AS A SINGLE DATASET
	ds = read_ar5_mon( path, variable, fn_prefix_filter, level=atmos_level )
	ds = read_ar5_mon( path, variable, filelist=files.tolist(), level=atmos_level, time_begin='2006-01-01', time_end='2100-12-31' )
	
	# generate climatology / anomalies
	climatology = ds.loc[ climatology_begin:climatology_end ].groupby( 'time.month' ).mean( 'time' )
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

	# LEGIT FILENAMING CONVENTION WITH A PREFIX / SUFFIX PATTERN REQUIRED HERE
	output_filenames = [ os.path.join( output_dir, 'hur_level11_akcan_' + str(i) + '.tif' ) for i in range( time_len ) ]

	# THIS ASSUMES THEY ARE THE ONLY FILES IN THE DIRECTORY -- COULD BE A GOTCHA
	cru_files = glob.glob( os.path.join( cru_path, '*.tif' ) )
	cru_files.sort()

	# this is effectively working
	cru_stack = [ rasterio.open( fn ).read( 1 ) for fn in cru_files ]
	cru_gen = cru_generator( 1271, cru_stack )

	# run in parallel
	pool = mp.Pool( 10 )
	out = pool.map( run, [{'src':src, 'output_filename':fn, 'dst':dst, 'cru':cru, 'src_crs':meta_4326[ 'crs' ], 'src_affine':meta_4326[ 'affine' ], \
							'dst_crs':meta_3338[ 'crs' ], 'dst_affine':meta_3338[ 'affine' ], 'dst_meta':meta_3338, 'operation':'add' } \
							for src,fn,cru in zip( np.vsplit( dat, time_len ), output_filenames, cru_gen ) ] ) 
	pool.close()



# this should be a launcher script of the downscaling and run across all of the data
# then we will run the conversions to vapor pressure on the outputs

# directory traversal <<- this stuff is now working 
import fnmatch
import functools
import itertools
import os
import pandas as pd
import xray, glob

def group_input_filenames( prefix, root_dir ):
	''' function that wraps some ugliness regarding returning the files we want to process '''
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

		for root_dir, dir_names, file_names in os.walk( path ):
			filter_partial = functools.partial(fnmatch.filter, file_names)

			for file_name in itertools.chain( *map( filter_partial, path_patterns ) ):
				yield os.path.join( root_dir, file_name )
	def version_grouper( x ):
		''' groupby function for grouping by filenames '''
		dir_path = os.path.dirname( x )
		fn, _ = os.path.splitext( os.path.basename( x ) )
		# remove dates from filename -- they have a hyphen
		fn_base = '_'.join([ i for i in fn.split( '_' ) if '-' not in i ])
		# return the path element that startswith 'v' this is the version attribute
		version = [ x for x in dir_path.split( os.path.sep ) if x.startswith( 'v' ) ]
		return '_'.join([ fn_base, version[0] ])
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

	# get all matches with prefix
	matches = pd.Series([ match for match in find_files( root_dir, [ prefix ] ) ])

	# group by version
	grouped = dict([ group for group in matches.groupby( matches.apply( version_grouper ) )])

	# group keys to DataFrame
	keys_df = pd.DataFrame({ key:key.split( '_' ) for key in grouped.keys() }).T

	# parse the keys / values and keep only latest versions
	keys_df_grouped = pd.concat([ drop_old_versions(i[1]) for i in keys_df.groupby( 2 ) ])

	# make a new dictionary holding the filenames grouped the way we want
	final_out = { k:v for k,v in grouped.iteritems() if k in keys_df_grouped.index.tolist() }
	return final_out


# lets try with models
base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data'
variables = [ 'tas', 'hur' ]
models = [ 'GISS-E2-R', 'IPSL-CM5A-LR' ] # fill this

# testing stuff
model = 'GISS-E2-R'
variable = 'hur'
prefix = variable + '_*' + model + '*'

for variable in variables:
	for model in models:
		prefix = variable + '_*' + model + '*'
		# run the input lineup
		tmp = group_input_filenames( prefix, base_path )
# run the above code here
os.system( 'python ' )


# # # # THIS IS A TESTING AREA 
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

# def open_dataset( fn ):

# def read_netcdfs(files, dim, transform_func=None):
#     def process_one_path(path):
#         # use a context manager, to ensure the file gets closed after use
#         with xray.open_dataset(path) as ds:
#             # transform_func should do some sort of selection or
#             # aggregation
#             if transform_func is not None:
#                 ds = transform_func(ds)
#             # load all data from the transformed dataset, to ensure we can
#             # use it after closing each original file
#             ds.load()
#             return ds

#     paths = sorted(glob(files))
#     datasets = [process_one_path(p) for p in paths]
#     combined = xray.concat(datasets, dim)
#     return combined


# def f( ds, time_begin, time_end ):
# 	''' transform function '''
# 	return ds[ 'hur' ].loc[ time_begin:time_end ][ :, 16, ... ]
 
os.path.join( path, fn_prefix_filter )

# model = 'GISS-E2-R'
variable = 'hur'
fn_prefix_filter = variable + '_*' + model + '*'
base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data'
files = group_input_filenames( fn_prefix_filter, base_path )
# files = files[ 'hur_Amon_rcp85_r1i1p1_v20121016' ]
# path = os.path.dirname(files.tolist()[0])
# paths = sorted(glob.glob(os.path.join( path, fn_prefix_filter )))

# parse the paths list to only have the data up to time_end

# concatenate
# tmp = xray.concat([ xray.open_dataset( i ).load() for i in paths[:5] ], 'time' )
# tmp.to_netcdf( 'test_output_xray.nc', mode='w', format='NETCDF4' )





