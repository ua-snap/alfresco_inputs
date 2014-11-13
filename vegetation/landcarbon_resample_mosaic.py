# resample the data to 1km and mosaick with the larger IEM extent
def reclassify( rasterio_rst, reclass_list, output_filename, band=1, creation_options=dict() ):
	'''
	MODIFIED: removed window walking...  too slow..

	this function will take a raster image as input and
	reclassify its values given in the reclass_list.
	The reclass list is a simple list of lists with the 
	following formatting:
		[[begin_range, end_range, new_value]]
		ie. [ [ 1,3,5 ],[ 3,4,6 ] ]
			* which converts values 1 to 2.99999999 to 5
				and values 3 to 3.99999999 to 6
				all other values stay the same.
	arguments:
		rasterio_rst = raster image instance from rasterio package
		reclass_list = list of reclassification values * see explanation
		band = integer marking which band you wnat to return from the raster
				default is 1.
		creation_options = gdal style creation options, but in the rasterio implementation
			* options must be in a dict where the key is the name of the gdal -co and the 
			  value is the value passed to that flag.  
			  i.e. 
			  	["COMPRESS=LZW"] becomes dict([('compress','lzw')])
	'''
	# this will update the metadata if a creation_options dict is passed as an arg.
	meta = rasterio_rst.meta
	if len( creation_options ) < 0:
		meta.update( creation_options )

	with rasterio.open( output_filename, mode='w', **meta ) as out_rst:
		band_arr = rasterio_rst.read_band( band ).data # this is a gotcha with the .data stuff
		for rcl in reclass_list:
			band_arr[ np.logical_and( band_arr >= rcl[0], band_arr < rcl[1] ) ] = rcl[2]
		out_rst.write_band( band, band_arr )
	return rasterio.open( output_filename )
def world2Pixel( geotransform, x, y ):
	"""
	Uses a geotransform (gdal.GetGeoTransform(), or rasterio equivalent)
	to calculate the pixel location of a geospatial coordinate
	"""
	ulX = geotransform[0]
	ulY = geotransform[3]
	xDist = geotransform[1]
	yDist = geotransform[5]
	rtnX = geotransform[2]
	rtnY = geotransform[4]
	pixel = int((x - ulX) / xDist)
	line = int((ulY - y) / xDist)
	return ( pixel, line )
def bounds_to_window( geotransform, rasterio_bounds ):
	'''
	return a rasterio window tuple-of-tuples used to read a subset
	of a rasterio raster file into a numpy array.  This is done by 
	passing the window argument in the:
		 dataset.read_band() or dataset.write_band()

	This function returns an object acceptable for use as a window 
	passed to the window argument.

	Notes:
	A window is a view onto a rectangular subset of a raster dataset
	and is described in rasterio by a pair of range tuples.
	window = ((row_start, row_stop), (col_start, col_stop))

	arguments:
		geotransform = 6-element rasterio transform 
			* typically from dataset.transform
		rasterio_bounds = (lower left x, lower left y, upper right x, upper right y)
			* typically from dataset.bounds in rasterio
	** This also requires the world2Pixel function.

	Depends:
		rasterio

	'''
	ll = rasterio_bounds[:2]
	ur = rasterio_bounds[2:]
	ll_xy, ur_xy = [ world2Pixel( geotransform, x, y ) for x, y in [ll, ur] ]
	return (( ur_xy[1], ll_xy[1]), ( ll_xy[0], ur_xy[0]))


if __name__ == '__main__':
	import os, rasterio, fiona
	import numpy as np

	input_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Output_Data'
	input_filename = os.path.join( input_dir, 'landcarbon_vegetation_modelinput_maritime_2001_v0_4.tif' )
	output_resampled = os.path.join( input_dir, 'landcarbon_vegetation_modelinput_maritime_2001_v0_4_1km.tif' )

	# if os.path.exists( output_resampled ):
	# 	os.remove( output_resampled )

	output_resampled = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/alaska_canada/alfresco_vegetation_mask_landcarbon.tif'

	# use gdalwarp to reproject 
	# command = 'gdalwarp -tr 1000 1000 -r mode -srcnodata None -dstnodata None -multi -co "COMPRESS=LZW" '+ input_filename + ' ' + output_resampled
	# os.system( command )

	maritime = rasterio.open( output_resampled )
	alfresco = rasterio.open( os.path.join( input_dir, 'alfresco_model_vegetation_input_2005.tif' ) )

	# reclassify martime to a common classification
	output_filename = output_resampled.replace( '.tif', '_iem_rcl.tif' )
	reclass_list = [[11, 12, 12], [10, 11, 11], [9, 10, 10], [8, 9, 9], [7, 8, 8], [1, 2, 0]]
	# [[1,2,0],[7,8,8],[8,9,9],[9,10,10],[10,11,11],[11,12,12]]
	maritime_rcl = reclassify( maritime, reclass_list, output_filename, band=1, creation_options={'compress':'lzw'} )

	# reclassify alfresco veg to a common classification
	output_filename = alfresco.name.replace( '.tif', '_iem_rcl.tif' )
	reclass_list = [[8,9,13],[9,10,8]]
	alfresco_rcl = reclassify( alfresco, reclass_list, output_filename, band=1, creation_options={'compress':'lzw'} )
	alfresco_rcl.close()

	# block out kodiak here.
	kodiak_mask_nalcms = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/shapefile_extents/kodiak_aoi_mask_nalcms.tif' )
	kodiak_shape_mask_nalcms = fiona.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/shapefile_extents/kodiak_aoi_mask_nalcms.shp' )
	output_filename = os.path.join( input_dir, 'mosaic_step1_mask_kodiak.tif' )
	'gdal_merge.py -o ' + output_filename + '-of GTiff -co "COMPRESS=LZW" -tap -v -init 255 -n 255 ' + alfresco_rcl.name + ' ' + kodiak_mask_nalcms.name

	alfresco_rcl = rasterio.open( output_filename ) # read it back in after masking

	# extend the extent of the martime data to match the extent of the alfresco data
	# overlay 2 maps and fill-in where North Pacific Maritime (class 13)
	window = bounds_to_window( alfresco_rcl.transform, maritime_rcl.bounds )
	meta = alfresco_rcl.meta
	meta.update( compress='lzw', nodata=255 )
	output_filename = os.path.join( input_dir, 'iem_model_vegetation_input_merged.tif' )
	with rasterio.open( output_filename, 'w', **meta ) as out:
		maritime_arr = maritime_rcl.read_band( 1 )
		akcan_arr = alfresco_rcl.read_band( 1, window=window )
		# akcan_arr[ (akcan_arr == 13) & (maritime_arr != 255) ] = maritime_arr
		np.place( akcan_arr, akcan_arr == 13, maritime_arr.ravel() )
		out.write_band( 1, akcan_arr, window=window )
		# extract kodiak from the maritime to pass into the new map
		window_kodiak = bounds_to_window( alfresco_rcl.transform, kodiak_shape_mask_nalcms.bounds )
		kodiak_arr = maritime_rcl.read_band( 1, window=window_kodiak )
		out.write_band( 1, kodiak_arr, window=window )







	# use polygon to convert Canada North Pacific Maritime to Upland Fores

	# use a polygon to set all other North Pacific Maritime to Out-of-bounds

	# pass in a colortable and write out




# [ TEM ] class changes from original

# All cells in the LandCover_iem_ALFRESCO_2005.tif file that were classed as 0 - No veg, were reset back to the original NALCMS values. 
# In this new file, cells classified as Cropland, Urban and Built-up, Water, and Snow and Ice were reclassified as 0 - No veg, 
# cells classified as Wetland were set to Wetland Tundra, and cells classified as Barren Lands were reclassified to Heath.

# no_veg = [15, 17, 18, 19]
# wetland_tundra = [14] 
# heath = [16]
