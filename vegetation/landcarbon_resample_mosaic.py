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
def hex_to_rgb( hex ):
	'''
	borrowed and modified from Matthew Kramer's blog:
		http://codingsimplicity.com/2012/08/08/python-hex-code-to-rgb-value/

	function to take a hex value and convert it into an RGB(A) representation.

	This is useful for generating color tables for a rasterio GTiff from a QGIS 
	style file (qml).  Currently tested for the QGIS 2.0+ style version.

	arguments:
		hex = hex code as a string

	returns:
		a tuple of (r,g,b,a), where the alpha (a) is ALWAYS 1.  This may need
		additional work in the future, but is good for the current purpose.
		** we need to figure out how to calculate that alpha value correctly.

	'''
	hex = hex.lstrip('#')
	hlen = len(hex)
	rgb = [ int( hex[i:i+hlen/3], 16 ) for i in range(0, hlen, hlen/3) ]
	rgb.insert(len(rgb)+1, 1)
	return rgb
def qml_to_ctable( qml ):
	'''
	take a QGIS style file (.qml) and converts it into a 
	rasterio-style GTiff color table for passing into a file.

	arguments:
		qml = path to a QGIS style file with .qml extension
	returns:
		dict of id as key and rgba as the values

	'''
	import xml.etree.cElementTree as ET
	tree = ET.ElementTree( file=qml  )
	return { int( i.get( 'value' ) ) : tuple( hex_to_rgb( i.get('color') ) ) for i in tree.iter( tag='item' ) }

if __name__ == '__main__':
	import os, rasterio, fiona, shutil, glob
	import numpy as np

	input_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Output_Data'
	input_filename = os.path.join( input_dir, 'landcarbon_vegetation_modelinput_maritime_2001_v0_4.tif' )
	output_resampled = os.path.join( input_dir, 'landcarbon_vegetation_modelinput_maritime_2001_v0_4_1km.tif' )

		# if os.path.exists( output_resampled ):
		# 	[ os.remove( i ) for i in glob.glob( output_resampled[:-3] + '*' ) ]

		# shutil.copy( os.path.join( input_dir, 'alfresco_model_vegetation_input_2005.tif' ), output_resampled )
		# maritime = rasterio.open( output_resampled ) # something odd here but this hack works -- does nothing
		# with rasterio.open( output_resampled, 'r+' ) as maritime:
		# 	arr = maritime.read_band( 1 ).data
		# 	arr[:] = 0
		# 	maritime.write_band( 1, arr )
		# 	del arr

		# # use gdalwarp to reproject and resample to the full akcanada domain
		# command = 'gdalwarp -r mode -multi -srcnodata None ' + input_filename + ' ' + output_resampled
		# os.system( command )



	# # now do the same thing with the combined maritime mask
	# output_resampled = os.path.join( input_dir, 'landcarbon_vegetation_modelinput_maritime_mask_1km.tif' )
	# shutil.copy( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/alaska_canada/alfresco_vegetation_mask.tif' , output_resampled )
	# maritime_mask = rasterio.open( output_resampled ) # something odd here but this hack works -- does nothing
	# with rasterio.open( output_resampled, 'r+' ) as maritime_mask:
	# 	arr = maritime_mask.read_band( 1 ).data
	# 	arr[:] = 0
	# 	maritime_mask.write_band( 1, arr )
	# 	del arr

	# combined_mask = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/combined_mask.tif'
	# command = 'gdalwarp -r mode  -multi -srcnodata None ' + combined_mask + ' ' + output_resampled
	# os.system( command )

	maritime = rasterio.open( output_resampled )
	alfresco = rasterio.open( os.path.join( input_dir, 'alfresco_model_vegetation_input_2005.tif' ) )

	# reclassify martime to a common classification
	output_filename = output_resampled.replace( '.tif', '_iem_rcl.tif' )
	reclass_list = [[11, 12, 12], [10, 11, 11], [9, 10, 10], [8, 9, 9], [7, 8, 8], [1, 2, 0]]
	maritime_rcl = reclassify( maritime, reclass_list, output_filename, band=1, creation_options={'compress':'lzw'} )
	maritime_rcl.close()

	# reclassify alfresco veg to a common classification
	output_filename = alfresco.name.replace( '.tif', '_iem_rcl.tif' )
	reclass_list = [[8,9,13],[9,10,8]]
	alfresco_rcl = reclassify( alfresco, reclass_list, output_filename, band=1, creation_options={'compress':'lzw'} )
	alfresco_rcl.close() # close it since there is no flush in rasterio
	alfresco_rcl = rasterio.open( output_filename ) # and reopen

	# alfresco_rcl = rasterio.open( output_filename ) # read it back in after masking
	maritime_rcl = rasterio.open( maritime_rcl.name )
	
	# generate an output colortable to pass to the new raster
	qml = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Output_Data/qgis_styles/iem_vegetation_modelinput_qgis_style.qml'
	cmap = qml_to_ctable( qml )

	maritime_mask_1k = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/combined_mask_alaska_only_akcan_1km.tif' )
	nlcd_saltwater_mask = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/nlcd_2001_land_cover_maritime_saltwater_mask_akcan_1km.tif' )
	maritime_mask_canada = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/combined_mask_canada_only_akcan_1km.tif' )
	temperate_rf_mask = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/digitized_removal_mask_temperate_rainforest_maritime.tif' )
	iem_mask = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/iem_mask_1km_akcan.tif' )

	# extend the extent of the martime data to match the extent of the alfresco data
	# overlay 2 maps and fill-in where North Pacific Maritime (class 13)
	meta = alfresco_rcl.meta
	meta.update( compress='lzw', nodata=255 )
	output_filename = os.path.join( input_dir, 'iem_model_vegetation_input_merged.tif' )
	with rasterio.open( output_filename, 'w', **meta ) as out:
		maritime_arr = maritime_rcl.read_band( 1 ).data
		maritime_mask = maritime_mask_1k.read_band(1) # ( maritime_arr.data > 0 ).astype( np.uint8 )
		maritime_mask_can_arr = maritime_mask_canada.read_band( 1 )
		akcan_arr = alfresco_rcl.read_band( 1 ).data
		nlcd_sw_arr = nlcd_saltwater_mask.read_band( 1 )
		temperate_rf_mask_arr = temperate_rf_mask.read_band( 1 )
		iem_mask_arr = iem_mask.read_band( 1 )

		# pass in values over the collective domain
		ind = np.where( (maritime_mask == 1) & (akcan_arr == 13) ) 
		akcan_arr[ ind ] = maritime_arr[ ind ]

		# remove the Canadian IEM domain Temperate Rainforest and convert to Maritime Upland Forest
		akcan_arr[ (maritime_mask_can_arr == 1) & (akcan_arr == 13) ] = 9
		# remove the rest of the errant Temperate Rainforest and convert to no veg
		akcan_arr[ (nlcd_sw_arr == 1) & (akcan_arr == 13) ] = 0
		akcan_arr[ (temperate_rf_mask_arr == 1) & (akcan_arr == 13) ] = 0
		
		# anything left on the canada side convert to Maritime Upland Forest
		akcan_arr[ (iem_mask_arr == 1) & (akcan_arr == 13) ] = 9

		# [not yet implemented] potentially find pixels in akcan that had veg data but do not in the new version
		#  and replace their value with their neighbors that are not nodata.
		
		out.write_band( 1, akcan_arr )
		out.write_colormap( 1, cmap )



# [ TEM ] class changes from original

# All cells in the LandCover_iem_ALFRESCO_2005.tif file that were classed as 0 - No veg, were reset back to the original NALCMS values. 
# In this new file, cells classified as Cropland, Urban and Built-up, Water, and Snow and Ice were reclassified as 0 - No veg, 
# cells classified as Wetland were set to Wetland Tundra, and cells classified as Barren Lands were reclassified to Heath.

# no_veg = [15, 17, 18, 19]
# wetland_tundra = [14] 
# heath = [16]

