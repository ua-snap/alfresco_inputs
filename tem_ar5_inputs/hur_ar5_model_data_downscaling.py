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


	# some setup pathing <<-- THIS TO BE CONVERTED TO ARGUMENTS AT COMMAND LINE
	modeled_fn = None
	historical_fn = None
	output_dir = '/home/UA/malindgren/Documents/hur/akcan/new'
	variable = 'hur'
	time_begin = '1850-01-01' # will change for future and historical
	time_end = '2005-12-31' # will change for future and historical
	climatology_begin = '1961-01-01'
	climatology_end = '1990-12-31'
	plev = 16 # this is also a None default!
	cru_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts20/akcan'

	# upack argparse'd inputs from the command line

	# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	# hardwired raster metadata meeting the ALFRESCO Model's needs for 
	# perfectly aligned inputs this is used as template metadata that 
	# is used in output generation. template raster filename below:
	# '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/
	#	TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
	meta_3338 = {'affine': Affine(2000.0, 0.0, -2173223.206087799, 
					0.0, -2000.0, 2548412.932644147),
				'count': 1,
				'crs': {'init':'epsg:3338'},
				'driver': u'GTiff',
				'dtype': 'float32',
				'height': 1186,
				'nodata': -3.4e+38,
				'transform': (-2173223.206087799, 2000.0, 0.0, 
					2548412.932644147,0.0,-2000.0),
				'width': 3218,
				'compress':'lzw'}
	# output template numpy array same dimensions as the template
	dst = np.empty( (1186, 3218) )
	
	# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	# WIRE THIS UP!
	# condition to deal with reading in historical data if needed.
	if modeled_fn is not None and historical_fn is not None:
		# read in both 
	elif historical_fn is not None:
		# read in historical
	else:
		NameError( 'ERROR: must have both modeled_fn and historical_fn, or just historical_fn' )

	# read in the prepped data:
	ds = xray.open_dataset( fn, chunks={'time':40} )

	# if there is a pressure level to extract, extract it
	if plev is not None:
		ds = ds[ :, plev, ... ]

	# generate climatology / anomalies
	climatology = ds.loc[ climatology_begin:climatology_end ].groupby( 'time.month' ).mean( 'time' )
	# RIGHT NOW THIS IS NOT THE PROPORTIONAL ANOMALY !!! THIS NEEDS TO BE A SWITCH TO DEAL WITH PROPORTION OR ABSOLUTE DIFF
	anomalies = ds.groupby( 'time.month' ) - climatology # anomalies = ds.groupby( 'time.month' ) / climatology

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


	# !we need some way to get some metadata from the NetCDF files we have read in...  or we need to assume things about the naming convention


	# LEGIT FILENAMING CONVENTION WITH A PREFIX / SUFFIX PATTERN REQUIRED HERE
	output_filenames = [ os.path.join( output_dir, '_'.join([variable, metric, model, scenario, experiment, month, year]) + '.tif' ) for i in range( time_len ) ]

	# THIS ASSUMES THEY ARE THE ONLY FILES IN THE DIRECTORY -- COULD BE A GOTCHA
	cru_files = glob.glob( os.path.join( cru_path, '*.tif' ) )
	cru_files.sort()

	# this is effectively working
	cru_stack = [ rasterio.open( fn ).read( 1 ) for fn in cru_files ]
	cru_gen = cru_generator( 1271, cru_stack ) # should 1271 be len(cru_files?)

	# run in parallel
	pool = mp.Pool( 10 )
	out = pool.map( run, [{'src':src, 'output_filename':fn, 'dst':dst, 'cru':cru, 'src_crs':meta_4326[ 'crs' ], 'src_affine':meta_4326[ 'affine' ], \
							'dst_crs':meta_3338[ 'crs' ], 'dst_affine':meta_3338[ 'affine' ], 'dst_meta':meta_3338, 'operation':'add' } \
							for src,fn,cru in zip( np.vsplit( dat, time_len ), output_filenames, cru_gen ) ] ) 
	pool.close()


