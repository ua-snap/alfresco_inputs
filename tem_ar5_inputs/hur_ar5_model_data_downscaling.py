# # # # #
# Tool to downscale the CMIP5 data from the PCMDI group. 
# # # # #

def cru_generator( n, cru_clim_list ):
	'''
	generator that will produce the cru climatologies with a
	generator and replicate for the total number of years in n
	'''
	months = [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12' ]
	for i in range( n ):
		for count, j in enumerate( cru_clim_list ):
			yield j
def standardized_fn_to_vars( fn ):
	''' take a filename string following the convention for this downscaling and break into parts and return a dict'''
	name_convention = [ 'variable', 'cmor_table', 'model', 'scenario', 'experiment', 'begin_time', 'end_time' ]
	fn = os.path.basename( fn )
	fn_list = fn.split( '.' )[0].split( '_' )
	return { i:j for i,j in zip( name_convention, fn_list )}
def downscale( src, dst, cru, src_crs, src_affine, dst_crs, dst_affine, output_filename, dst_meta, \
		method='cubic_spline', operation='add', output_dtype='float32', **kwargs ):
	'''
	operation can be one of two keywords for the operation to perform the delta downscaling
	- keyword strings are one of: 'add'= addition, 'mult'=multiplication, or 'div'=division (not implemented)
	- method can be one of 'cubic_spline', 'nearest', 'bilinear' and must be input as a string.
	- output_dtype can be one of 'int32', 'float32'
	'''
	from rasterio.warp import reproject, RESAMPLING
	def add( cru, anom ):
		return cru + anom
	def mult( cru, anom ):
		return cru * anom
	def div( cru, anom ):
		# return cru / anom
		# this one may not be useful, but the placeholder is here 
		return NotImplementedError

	# switch to deal with numeric output dtypes
	dtypes_switch = {'int32':np.int32, 'float32':np.float32}

	# switch to deal with different resampling types
	method_switch = { 'nearest':RESAMPLING.nearest, 'bilinear':RESAMPLING.bilinear, 'cubic_spline':RESAMPLING.cubic_spline }
	method = method_switch[ method ]

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
		out.write_band( 1, downscaled.astype( dtypes_switch[ output_dtype ] ) )
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

	# parse the commandline arguments
	parser = argparse.ArgumentParser( description='preprocess cmip5 input netcdf files to a common type and single files' )
	parser.add_argument( "-mi", "--modeled_fn", nargs='?', const=None, action='store', dest='modeled_fn', type=str, help="path to modeled input filename (NetCDF); default:None" )
	parser.add_argument( "-hi", "--historical_fn", nargs='?', const=None, action='store', dest='historical_fn', type=str, help="path to historical input filename (NetCDF); default:None" )

	parser.add_argument( "-o", "--output_dir", action='store', dest='output_dir', type=str, help="string path to the output folder containing the new downscaled outputs" )
	parser.add_argument( "-bt", "--begin_time", action='store', dest='begin_time', type=str, help="string in format YYYYMM of the beginning month/year" )
	parser.add_argument( "-et", "--end_time", action='store', dest='end_time', type=str, help="string in format YYYYMM of the ending month/year" )
	parser.add_argument( "-cbt", "--climatology_begin_time", nargs='?', const='196101', action='store', dest='climatology_begin', type=str, help="string in format YYYYMM or YYYY of the beginning month and potentially (year) of the climatology period" )
	parser.add_argument( "-cet", "--climatology_end_time", nargs='?', const='199012', action='store', dest='climatology_end', type=str, help="string in format YYYYMM or YYYY of the ending month and potentially (year) of the climatology period" )
	parser.add_argument( "-plev", "--plev", nargs='?', const=None, action='store', dest='plev', type=int, help="integer value (in millibars) of the desired pressure level to extract, if there is one." )
	parser.add_argument( "-cru", "--cru_path", action='store', dest='cru_path', type=str, help="path to the directory storing the cru climatology data derived from CL2.0" )
	parser.add_argument( "-at", "--anomalies_calc_type", nargs='?', const='absolute', action='store', dest='anomalies_calc_type', type=int, help="string of 'proportional' or 'absolute' to inform of anomalies calculation type to perform." )

	# parse and unpack args
	args = parser.parse_args()

	# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	# [NOTE]: hardwired raster metadata meeting the ALFRESCO Model's needs for 
	# perfectly aligned inputs this is used as template metadata that 
	# is used in output generation. template raster filename below:
	# '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/
	#	TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
	meta_3338 = {'affine': A(2000.0, 0.0, -2173223.206087799, 
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
	# condition to deal with reading in historical data if needed.
	if modeled_fn is not None and historical_fn is not None:
		# read in both
		ds = xray.open_dataset( modeled_fn )
		ds = ds[ variable ].load()
		clim_ds = xray.open_dataset( historical_fn )
		clim_ds = clim_ds[ variable ].load()
		# generate climatology / anomalies
		clim_ds = clim_ds.loc[ {'time':slice(climatology_begin,climatology_end)} ]
		climatology = clim_ds.groupby( 'time.month' ).mean( 'time' )
		# parse the input name for some file metadata
		output_naming_dict = standardized_fn_to_vars( modeled_fn )
	elif historical_fn is not None:
		# read in historical
		ds = xray.open_dataset( historical_fn )
		ds = ds[ variable ].load()
		# generate climatology / anomalies
		climatology = ds.loc[ {'time':slice(climatology_begin,climatology_end)} ]
		climatology = climatology.groupby( 'time.month' ).mean( 'time' )
		# parse the input name for some file metadata
		output_naming_dict = standardized_fn_to_vars( historical_fn )
	else:
		NameError( 'ERROR: must have both modeled_fn and historical_fn, or just historical_fn' )

	# if there is a pressure level to extract, extract it
	if plev is not None:
		plevel, = np.where( ds.plev == plev )
		ds = ds[ :, plevel[0], ... ]
		climatology = climatology[ :, plevel[0], ... ]

	# deal with different anomaly calculation types
	if anomalies_calc_type == 'absolute':
		anomalies = ds.groupby( 'time.month' ) - climatology
	elif anomalies_calc_type == 'proportional':
		anomalies = ds.groupby( 'time.month' ) / climatology
	else:
		NameError( 'anomalies_calc_type can only be one of "absolute" or "proportional"' )

	# some setup of the output raster metadata
	time_len, rows, cols = anomalies.shape
	crs = 'epsg:4326'
	affine = A( *[np.diff( ds.lon )[ 0 ], 0.0, -180.0, 0.0, -np.diff( ds.lat )[ 0 ], 90.0] )
	count = time_len
	resolution = ( np.diff( ds.lat )[ 0 ], np.diff( ds.lon )[ 0 ] )

	# shift the grid to Greenwich Centering
	dat, lons = shiftgrid( 180., anomalies[:], anomalies.lon.data, start=False )

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

	# build some filenames for the outputs to be generated
	months = [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12' ]
	years = [ str(year) for year in range( int(time_begin[:4]), int(time_end[:4]) + 1, 1 ) ]
	# combine the months and the years
	combinations = [ (month, year) for year in years for month in months ]
	output_filenames = [ os.path.join( output_dir, '_'.join([output_naming_dict['variable'], 'metric', output_naming_dict['model'], output_naming_dict['scenario'], output_naming_dict['experiment'], month, year]) + '.tif' ) for month, year in combinations ]

	# load the baseline CRU CL2.0 data 
	# [NOTE]: THIS ASSUMES THEY ARE THE ONLY FILES IN THE DIRECTORY -- COULD BE A GOTCHA
	cru_files = glob.glob( os.path.join( cru_path, '*.tif' ) )
	cru_files.sort()
	cru_stack = [ rasterio.open( fn ).read( 1 ) for fn in cru_files ]
	cru_gen = cru_generator( len(output_filenames), cru_stack )

	# cleanup some uneeded vars that are hogging RAM
	del clim_ds, climatology, ds, anomalies

	# run in parallel using PATHOS
	pool = mp.Pool( 8 )
	out = pool.map( run, [{'src':src, 'output_filename':fn, 'dst':dst, 'cru':cru, 'src_crs':meta_4326[ 'crs' ], 'src_affine':meta_4326[ 'affine' ], \
							'dst_crs':meta_3338[ 'crs' ], 'dst_affine':meta_3338[ 'affine' ], 'dst_meta':meta_3338, 'operation':'mult' } \
							for src,fn,cru in zip( np.vsplit( dat, time_len ), output_filenames, cru_gen ) ] ) 
	pool.close()


# # # # # # # # #SOME TESTING AND EXAMPLE GENERATION AREA # # # # # # # # # 
# # some setup pathing <<-- THIS TO BE CONVERTED TO ARGUMENTS AT COMMAND LINE
# historical_fn = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped/GFDL-CM3/hur/hur_Amon_GFDL-CM3_historical_r1i1p1_186001_200512.nc' 
# modeled_fn = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped/GFDL-CM3/hur/hur_Amon_GFDL-CM3_rcp26_r1i1p1_200601_210012.nc' 

# variable = 'hur'
# metric = 'pct'
# output_dir = '/home/UA/malindgren/Documents/hur/akcan/new'
# time_begin = '2006-01' # will change for future and historical
# time_end = '2100-12' # will change for future and historical
# climatology_begin = '1961'
# climatology_end = '1990'
# plev = 1000 # this is in millibar data, this is also a None default!
# cru_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts20/akcan'
# anomalies_calc_type = 'proportional'


