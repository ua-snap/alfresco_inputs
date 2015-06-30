
import numpy as np # hack to solve a lib issue in the function args of xyztogrid

def cru_xyz_to_shp( in_xyz, lon_col, lat_col, crs, output_filename ):
	'''
	convert the cru cl2.0 1961-1990 Climatology data to a shapefile.
		*can handle the .dat format even if compressed with .gzip extension.

	PARAMETERS:
	-----------
	in_xyz = path to the .dat or .dat.gz downloaded cru cl2.0 file from UK Met Office site
	lon_col = string name of column storing longitudes
	lat_col = string name of column storing latitudes
	crs = proj4string or epsg code
	output_filename = string path to the output filename to be created

	RETURNS
	-------
	output_filename as string

	'''
	colnames = ['lat', 'lon', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
	from shapely.geometry import Point
	import pandas as pd
	import geopandas as gpd
	import os

	if os.path.splitext( in_xyz )[1] == '.gz':
		cru_df = pd.read_csv( in_xyz, delim_whitespace=True, compression='gzip', header=None, names=colnames )
	else:
		cru_df = pd.read_csv( in_xyz, delim_whitespace=True, header=None, names=colnames )

	# create a column named geometry with shapely geometry objects for each row
	def f( x ):
		''' return a Point shapely object for each x,y pair'''
		return Point( x.lon, x.lat )
	
	cru_df[ 'geometry' ] = cru_df.apply( f, axis=1 )
	cru_df = gpd.GeoDataFrame( cru_df ) # convert to GeoDataFrame
	cru_df.to_file( output_filename, 'ESRI Shapefile' )
	return output_filename
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
def crop_to_bounds( rasterio_rst, bounds, output_filename, mask=None, mask_value=None ):
	'''
	take a rasterio raster object and crop it to a smaller bounding box
	masking is supported where masked values are 0 and unmasked values are 1

	PARAMETERS
	----------
	rasterio_rst = rasterio raster object
	bounds = rasterio style bounds (left, bottom, right, top)
	output_filename = string path to the raster file to be created
	mask = a 2d numpy array of the same shape as rasterio_rst with
		masked values = 0 and unmasked = 1

	RETURNS
	-------
	file path to the newly created file -- essentially the value of output_filename

	'''
	from rasterio import Affine as A
	window = rasterio_rst.window( *bounds )
	xmin, ymin, xmax, ymax = rasterio_rst.window_bounds( window )
	row_res, col_res = rasterio_rst.res
	arr = rasterio_rst.read( 1, window=window )
	
	if mask:
		arr[ mask != 1 ] = mask_value
		nodata = mask_value
	else:
		nodata = rasterio_rst.meta[ 'nodata' ]

	meta = {}
	meta.update( compress='lzw',
				affine=A( col_res, 0.0, xmin, 0.0, -row_res, ymax ),
				height=row_res,
				width=col_res,
				transform=[xmin, col_res, 0.0, ymax, 0.0, -row_res],
				crs=rasterio_rst.meta,
				nodata=nodata,
				dtype=rasterio_rst.meta[ 'dtype' ],
				count=1,
				driver=u'GTiff' )

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write_band( 1, arr )
	return output_filename

if __name__ == '__main__':
	import os, rasterio, glob, fiona
	import numpy as np
	import pandas as pd
	import geopandas as gpd
	from rasterio import Affine as A

	base_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data'

	# open the Climatic Research Unit (CRU) CL2.0 data downloaded from:
	#		 http://www.cru.uea.ac.uk/cru/data/hrg/tmc/
	cru_filename = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS20/grid_10min_reh.dat.gz'
	cru_path = os.path.join( base_path, 'cru_ts20' )

	# read in the gzipped .dat file downloaded from the MET Office UK
	colnames = [ 'lat', 'lon', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12' ]
	cru_df = pd.read_csv( cru_filename, delim_whitespace=True, compression='gzip', header=None, names=colnames )

	# convert to point shapefile
	cru_shp_fn = os.path.join( cru_path, 'cru_ts20_1961_1990_climatology.shp' )
	cru_xyz_to_shp( cru_filename, 'lon', 'lat', {'init':'epsg:4326'}, cru_shp_fn )
	
	# template dataset
	akcan_template = rasterio.open( os.path.join( base_path, 'templates', 'tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif' ) )
	resolution = akcan_template.res
	akcan_meta = akcan_template.meta

	# pad the bounds of the akcan template dataset
	crs = { 'init':'epsg:3338' }
	akcan_ext_fn = os.path.join( base_path, 'extents', 'akcan_extent.shp' )
	npixels = ( -200, -2000, 200, 200 )
	pad_bounds( akcan_template, npixels, crs, akcan_ext_fn )

	# filename for a newly clipped and reprojected shapefile using the above padded bounds shape
	expanded_ext_fn = os.path.join( base_path, 'intermediate', 'cru_ts20_1961_1990_climatology_3338_akcan_expanded.shp' )

	# reproject / crop to the AKCAN extent, the cru shapefile built above using ogr2ogr
	os.system( "ogr2ogr -progress -wrapdateline -overwrite -f 'ESRI Shapefile' -clipdst " + akcan_ext_fn +  " -s_srs 'EPSG:4326' -t_srs 'EPSG:3338' " + expanded_ext_fn + " " + cru_shp_fn )

	# generate metadata for the expanded extent to interpolate to
	xmin, ymin, xmax, ymax = fiona.open( akcan_ext_fn ).bounds
	cols = (xmax - xmin) / resolution[1]
	rows = (ymax - ymin) / resolution[0]

	# copy/update metadata to expanded extent
	expanded_akcan_meta = akcan_meta
	expanded_akcan_meta[ 'affine' ] = A( resolution[0], 0.0, xmin, 0.0, -resolution[1], ymax )
	expanded_akcan_meta[ 'crs' ] = { 'init':'epsg:3338' }
	expanded_akcan_meta[ 'height' ] = rows
	expanded_akcan_meta[ 'width' ] = cols
	expanded_akcan_meta[ 'transform' ] = expanded_akcan_meta[ 'affine' ].to_gdal()

	# read in the clipped and reprojected cru shapefile using geopandas
	cru_df_akcan = gpd.read_file( expanded_ext_fn )

	# update lon and lat to the 3338
	cru_df_akcan.lon = cru_df_akcan.geometry.apply( lambda x: x.x )
	cru_df_akcan.lat = cru_df_akcan.geometry.apply( lambda x: x.y )

	# build the interpolation input values
	x = np.array(cru_df_akcan.lon.tolist())
	y = np.array(cru_df_akcan.lat.tolist())

	# build the output grid
	xi = np.linspace( xmin, xmax, cols )
	yi = np.linspace( ymin, ymax, rows )
	xi, yi = np.meshgrid( xi, yi )

	# build some args	
	months = ['01','02','03','04','05','06','07','08','09','10','11','12']
	output_filenames = [ os.path.join( akcan_path, 'hur_cru_cl20_akcan_'+month+'_1961_1990.tif' ) for month in months ]

	def interpolate_akcan( x, y, z, grid, expanded_meta, template_rst, output_filename, mask=None, mask_value=None, method='cubic', output_dtype=np.float32 ):
		'''
		interpolate across the alaska canada domains and crop / mask to that extent
		'''	
		cru_interp = xyz_to_grid( x, y, z, grid, method='cubic', output_dtype=np.float32 )
		cru_interp = np.nan_to_num( cru_interp )

		# convert to in memory rasterio object
		expanded_meta.update( driver='MEM' )
		cru_interpolated = rasterio.open( '', mode='w', **expanded_meta )
		cru_interpolated.write_band( 1, cru_interp )

		akcan = crop_to_bounds( cru_interpolated, template_rst.bounds, output_filename, mask, mask_value )
		
		meta = template_rst.meta
		meta.update( compress='lzw' )
		with rasterio.open( output_filename, 'w', **meta ) as out:
			mask = template_rst.read_mask()
			akcan[ mask == 0 ] = meta[ 'nodata' ]
			out.write_band( 1, akcan )
		return output_filename

	def run( args ):
		return interpolate_akcan( **args )

	# run it in parallel
	from pathos import multiprocessing as mp
	args_list = [ { 'x':x, 'y':y, 'z':cru_df_akcan[ month ], 'grid':(xi,yi), 'expanded_meta':expanded_akcan_meta, 'template_rst':akcan_template, 'output_filename':out_fn } for month, out_fn in zip( months, output_filenames ) ]
	pool = mp.Pool( 10 )
	pool.map( run, args_list )
	pool.close()


	

	# # interpolate points to grid
	# cru_interpolated = [ xyz_to_grid( x, y, np.array(cru_df_akcan[ month ].tolist()), (xi, yi), method='cubic', output_dtype=np.float32 ) for month in months ]

	# # crop it to the akcan extent
	# [ crop_to_bounds( rasterio_rst, akcan_template.bounds, output_filename, mask=None, mask_value ) for rst in cru_interpolated ]
	

	# akcan_path = os.path.join( base_path, 'akcan' )
	# if not os.path.exists( akcan_path ):
	# 	os.makedirs( akcan_path )

	# # set up some output raster filenames


	
	# meta.update( compress='lzw' )
	# with rasterio.open( output_filename, 'w', **meta ) as out:
	# 	out.write_band( 1, zi )


