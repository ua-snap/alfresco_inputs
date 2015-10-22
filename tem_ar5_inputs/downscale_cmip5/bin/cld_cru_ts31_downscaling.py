import numpy as np

def shiftgrid(lon0,datain,lonsin,start=True,cyclic=360.0):
	import numpy as np
	"""
	Shift global lat/lon grid east or west.
	.. tabularcolumns:: |l|L|
	==============   ====================================================
	Arguments        Description
	==============   ====================================================
	lon0             starting longitude for shifted grid
					 (ending longitude if start=False). lon0 must be on
					 input grid (within the range of lonsin).
	datain           original data with longitude the right-most
					 dimension.
	lonsin           original longitudes.
	==============   ====================================================
	.. tabularcolumns:: |l|L|
	==============   ====================================================
	Keywords         Description
	==============   ====================================================
	start            if True, lon0 represents the starting longitude
					 of the new grid. if False, lon0 is the ending
					 longitude. Default True.
	cyclic           width of periodic domain (default 360)
	==============   ====================================================
	returns ``dataout,lonsout`` (data and longitudes on shifted grid).
	"""
	if np.fabs(lonsin[-1]-lonsin[0]-cyclic) > 1.e-4:
		# Use all data instead of raise ValueError, 'cyclic point not included'
		start_idx = 0
	else:
		# If cyclic, remove the duplicate point
		start_idx = 1
	if lon0 < lonsin[0] or lon0 > lonsin[-1]:
		raise ValueError('lon0 outside of range of lonsin')
	i0 = np.argmin(np.fabs(lonsin-lon0))
	i0_shift = len(lonsin)-i0
	if np.ma.isMA(datain):
		dataout  = np.ma.zeros(datain.shape,datain.dtype)
	else:
		dataout  = np.zeros(datain.shape,datain.dtype)
	if np.ma.isMA(lonsin):
		lonsout = np.ma.zeros(lonsin.shape,lonsin.dtype)
	else:
		lonsout = np.zeros(lonsin.shape,lonsin.dtype)
	if start:
		lonsout[0:i0_shift] = lonsin[i0:]
	else:
		lonsout[0:i0_shift] = lonsin[i0:]-cyclic
	dataout[...,0:i0_shift] = datain[...,i0:]
	if start:
		lonsout[i0_shift:] = lonsin[start_idx:i0+start_idx]+cyclic
	else:
		lonsout[i0_shift:] = lonsin[start_idx:i0+start_idx]
	dataout[...,i0_shift:] = datain[...,start_idx:i0+start_idx]
	return dataout,lonsout
def bounds_to_extent( bounds ):
	'''
	take input rasterio bounds object and return an extent
	'''
	l,b,r,t = bounds
	return [ (l,b), (r,b), (r,t), (l,t), (l,b) ]
def padded_bounds( rst, npixels, crs ):
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
	return new_bounds
	
def xyz_to_grid( x, y, z, grid, method='cubic', output_dtype=np.float32 ):
	'''
	interpolate points to a grid. simple wrapper around
	scipy.interpolate.griddata. Points and grid must be
	in the same coordinate system
	x = 1-D np.array of x coordinates / x,y,z must be same length
	y = 1-D np.array of y coordinates / x,y,z must be same length
	z = 1-D np.array of z coordinates / x,y,z must be same length
	grid = tuple of meshgrid as made using numpy.meshgrid()
			order (xi, yi)
	method = one of 'cubic', 'near', linear
	'''
	import numpy as np
	from scipy.interpolate import griddata

	zi = griddata( (x, y), z, grid, method=method )
	zi = np.flipud( zi.astype( output_dtype ) )
	return zi
def run( args ):
	''' 
	simple function wrapper for unpacking an argument dict 
	to the downscale function for getting around the single 
	argument pass to multiprocessing.map implementation issue.
	'''
	return( downscale( **args ) )

if __name__ == '__main__':
	import rasterio, xray, os, glob
	import geopandas as gpd
	import pandas as pd
	import numpy as np
	from collections import OrderedDict
	from shapely.geometry import Point
	from pathos import multiprocessing as mp

	# filenames
	base_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data'
	cld_ts31 = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS31/cru_ts_3_10.1901.2009.cld.dat.nc'
	template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
	# this is the set of modified GTiffs produced in the conversion procedure with the ts2.0 data
	cld_ts20 = '' # read in the already pre-produced files.  They should be in 10'...  or maybe I need to change that.
	climatology_begin = '1961'
	climatology_end = '1990'

	# open with xray
	cld_ts31 = xray.open_dataset( cld_ts31 )

	# calculate the anomalies
	clim_ds = cld_ts31.loc[ {'time':slice(climatology_begin,climatology_end)} ]
	climatology = clim_ds.cld.groupby( 'time.month' ).mean( 'time' )
	anomalies = cld_ts31.cld.groupby( 'time.month' ) / climatology

	# rotate the anomalies to pacific centered latlong -- this is already in the greenwich latlong
	dat, lons = shiftgrid( 0., anomalies, anomalies.lon )

	# # generate an expanded extent (from the template_raster) to interpolate across
	# template_raster = rasterio.open( template_raster_fn )
	# output_resolution = (1000.0, 1000.0) # hardwired, but we are building this for IEM which requires 1km
	# template_meta = template_raster.meta

	# npixels = ( -200, -2000, 200, 200 )
	# expanded_bounds = padded_bounds( template_raster, npixels, template_raster.crs )

	# # interpolate to a new grid
	# get longitudes and latitudes using meshgrid
	lo, la = [ i.ravel() for i in np.meshgrid( lons, cld_ts31.lat ) ]
	
	# convert into GeoDataFrame and drop all the NaNs
	df_list = [ pd.DataFrame({ 'anom':i.ravel(), 'lat':la, 'lon':lo }).dropna( axis=0, how='any' ) for i in dat ]
	xi, yi = np.meshgrid( lons, cld_ts31.lat )
	
	if __name__ == '__main__':
		# interpolate to the 0.5 degree grid to fill in some areas along the always problematic coastline
		pool = mp.Pool( 32 )
		out = pool.map( lambda df: xyz_to_grid( np.array(df['lon'].tolist()), np.array(df['lat'].tolist()), np.array(df['anom'].tolist()), grid=(xi,yi), method='cubic' ) , df_list[:2] )
		pool.close()

	# then restack these things
	anomalies_interp = np.rollaxis( np.dstack( out ), -1 )

	# rotate back into the Greenwich Centered Latlong
	anomalies_greenwich, lons = shiftgrid( 180., anomalies_interp, lons, start=False )
	
	# for reasons I do not currently understand (or care to figure out) I need to shift these new lons by 0.5 degrees.
	lons = lons + 0.5

	# cleanup some big objects
	del out, df_list, anomalies_interp
	
	# downscale - which will involves reprojection to the akcan extent and multiplying the relative anomalies to the 
	#  baseline raster data from CRU CL2.0

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
		out = np.zeros( dst.shape ) 
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
			out.write_band( 1, downscaled.astype( dtypes_switch[ output_dtype ] ).data )
		return output_filename



downscale( src, dst, cru, src_crs, src_affine, dst_crs, dst_affine, output_filename, dst_meta, method='cubic_spline', operation='mult', output_dtype='float32', **kwargs )



	# make a new xray.Dataset with these outputs from the anomalies calculation for interpolatiobn
	# anomalies_flip = xray.Dataset( {'cld_anom': (('time', 'lat', 'lon'), anomalies_interp)}, {'time':cld_ts31.time, 'lon':lons, 'lat':cld_ts31.lat } )

	# maybe we should just interpolate it to the final extent we want from here?
	# then just loop through the series of the sunp 10' climatologies? or just interp right to the 1km?




# # # # # # # # # # #  TESTING SOME NEW CODE 	# # # # # # # # # # #
# pad the bounds of the akcan template dataset
# crs = { 'init':'epsg:3338' } # this needs to be gotten from the file itself
# extent_path = os.path.join( base_path, 'extents_TEST' )
# if not os.path.exists( extent_path ):
# 	os.makedirs( extent_path )
# new_ext_fn = os.path.join( extent_path, 'akcan_extent_TEST.shp' )
npixels = ( -200, -2000, 200, 200 )
expanded_bounds = padded_bounds( template_raster, npixels, template_raster.crs )

# we should use the above padded bounds to filter the points 





def expand_templateds_extent( template_raster_fn, expansion_tuple=(0,0,0,0) ):
	'''
	args:
		template_raster_fn = [ str ] path to the template raster to match to
		expansion_tuple = [ int ] 4 element tuple of (left, bottom, right, top) indicating how many pixels
			to expand in each direction of the input template raster.

	'''
# template dataset
template_raster = rasterio.open( template_raster_fn )
resolution = template_raster.res
template_meta = template_raster.meta

# pad the bounds of the akcan template dataset
# crs = { 'init':'epsg:3338' } # this needs to be gotten from the file itself
extent_path = os.path.join( cru_path, 'extents' )
if not os.path.exists( extent_path ):
	os.makedirs( extent_path )
new_ext_fn = os.path.join( extent_path, 'akcan_extent.shp' )
npixels = ( -200, -2000, 200, 200 )
pad_bounds( template_raster, npixels, template_raster.crs, new_ext_fn )

# filename for a newly clipped and reprojected shapefile using the above padded bounds shape
intermediate_path = os.path.join( cru_path, 'intermediate' )
if not os.path.exists( intermediate_path ):
	os.makedirs( intermediate_path )

expanded_ext_fn = os.path.join( intermediate_path, variable + '_cru_ts20_1961_1990_climatology_3338_akcan_expanded.shp' )

# reproject / crop to the AKCAN extent, the cru shapefile built above using ogr2ogr
os.system( "ogr2ogr -overwrite -f 'ESRI Shapefile' -clipdst " + new_ext_fn +  " -s_srs 'EPSG:4326' -t_srs 'EPSG:3338' " + \
			expanded_ext_fn + " " + cru_shp_fn )
# -wrapdateline -- removed since it is not a geog srs output

# generate metadata for the expanded extent to interpolate to
xmin, ymin, xmax, ymax = fiona.open( new_ext_fn ).bounds
cols = (xmax - xmin) / resolution[1]
rows = (ymax - ymin) / resolution[0]

# copy/update metadata to expanded extent
expanded_meta = template_meta
expanded_meta[ 'affine' ] = A( resolution[0], 0.0, xmin, 0.0, -resolution[1], ymax )
expanded_meta[ 'crs' ] = { 'init':'epsg:3338' }
expanded_meta[ 'height' ] = rows
expanded_meta[ 'width' ] = cols
expanded_meta[ 'transform' ] = expanded_meta[ 'affine' ].to_gdal()





