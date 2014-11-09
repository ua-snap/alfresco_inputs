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

if __name__ == '__main__':
	import os, rasterio, fiona
	import numpy as np

	input_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Output_Data'
	input_filename = os.path.join( input_dir, 'landcarbon_vegetation_modelinput_maritime_2001_v0_4.tif' )
	output_resampled = os.path.join( input_dir, 'landcarbon_vegetation_modelinput_maritime_2001_v0_4_1km.tif' )

	if os.path.exists( output_resampled ):
		os.remove( output_resampled )

	# use gdalwarp to reproject 
	command = 'gdalwarp -tr 1000 1000 -r mode -srcnodata None -dstnodata None -multi -co "COMPRESS=LZW" '+ input_filename + ' ' + output_resampled
	os.system( command )

	maritime = rasterio.open( output_resampled )
	alfresco = rasterio.open( os.path.join( input_dir, 'alfresco_model_vegetation_input_2005.tif' ) )

	# reclassify martime to a common classification
	output_filename = output_resampled.replace( '.tif', '_iem_rcl.tif' )
	reclass_list = [[11, 12, 12], [10, 11, 11], [9, 10, 10], [8, 9, 9], [7, 8, 8], [1, 2, 0]]
	[[1,2,0],[7,8,8],[8,9,9],[9,10,10],[10,11,11],[11,12,12]]
	maritime_rcl = reclassify( maritime, reclass_list, output_filename, band=1, creation_options={'compress':'lzw'} )

	# reclassify alfresco veg to a common classification
	output_filename = alfresco.name.replace( '.tif', '_iem_rcl.tif' )
	reclass_list = [[8,9,13],[9,10,8]]
	alfresco_rcl = reclassify( alfresco, reclass_list, output_filename, band=1, creation_options={'compress':'lzw'} )

	# extend the extent of the martime data to match the extent of the alfresco data

	# overlay 2 maps and fill-in where North Pacific Maritime (class 13)

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
