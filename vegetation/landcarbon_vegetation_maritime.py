# # # # 
# reclassify NLCD 2001 Landcover for the LandCarbon project 
# South-Central Alaska mountains domain
# # # # 

#  Legend: 
#  0 - NoVeg 
#  1 - Black Spruce 
#  2 - White Spruce 
#  3 - Deciduous 
#  4 - Shrub Tundra 
#  5 - Graminoid  Tundra 
#  6 - Wetland Tundra 
#  7 - Barren lichen-moss 
#  8 - Temperate Rainforest (there will still be some in Canada to Contend with)
#  9 - Heath
#  10 - Upland
#  11 - Forested Wetland
#  12 - Fen
#  13 - Other Veg
#  14 - Alpine (optional)
#  15 - Alder

#  ** 16 no veg
#  ** 17 saltwater

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
		band_arr = rasterio_rst.read_band( band )
		for rcl in reclass_list:
			band_arr[ np.logical_and( band_arr >= rcl[0], band_arr < rcl[1] ) ] = rcl[2]
		out_rst.write_band( band, band_arr )
	return rasterio.open( output_filename )


if __name__ == '__main__':
	import os, rasterio, fiona, shutil
	from rasterio import features
	from rasterio.warp import reproject, RESAMPLING
	import numpy as np

	# some initial setup
	version_num = 'v0_4'
	input_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/seak'
	output_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Output_Data'
	os.chdir( output_dir )
	meta_updater = dict( driver='GTiff', dtype=rasterio.int16, compress='lzw', crs={'init':'epsg:3338'}, count=1, nodata=None )
	
	input_paths = { 
			'lc01':os.path.join( input_dir, 'NLCD_land_cover_AKNPLCC.tif' ),
			'cp01':os.path.join( input_dir, 'NLCD_pct_canopy_AKNPLCC.tif' ),
			'logged':os.path.join( input_dir, 'AKNPLCC_2ndGrowth.tif' ),
			'seak_mask':os.path.join( input_dir, 'seak_aoi.tif' ),
			'scak_mask':os.path.join( input_dir, 'scak_aoi.tif' )
	}

	# # duplicate the input lc map to modify in-place
	output_filename = os.path.join( output_dir, 'landcarbon_model_vegetation_input_maritime_2001_' + version_num + '_tmp_duplicate.tif' )
	shutil.copy( input_paths[ 'lc01' ], output_filename )

	# mask the non-seak region
	with rasterio.open( output_filename, mode='r+' ) as lc:
		lc_mod = lc.read_band( 1 )
		lc_mod_mask = lc_mod.mask
		seak_mask = rasterio.open( input_paths[ 'seak_mask' ] ).read_band( 1 )
		# seak_mask = np.ma.masked_where( seak_mask==0, seak_mask )
		# cur_mask = np.ma.mask_or( seak_mask.mask, lc_mod_mask )
		# lc_mod = np.ma.masked_array( lc_mod.data, cur_mask )

		# collapse initial undesired classes to noveg
		lc_mod = lc.read_band( 1 )
		lc_mod[ (lc_mod >= 0) & (lc_mod <= 31 ) & (seak_mask == 1) ] = 1 # potentially 0 later on

		# upland forest / fen
		canopy = rasterio.open( input_paths[ 'cp01' ] ).read_band( 1 )
		lc_mod[ (lc_mod == 42) & (canopy > 20) & (seak_mask == 1) ] = 10
		lc_mod[ (lc_mod == 42) & (canopy <= 20) & (seak_mask == 1) ] = 12

		# alder / shrubland
		lc_mod[ (lc_mod == 81) | (lc_mod == 82) & (seak_mask == 1) ] = 15

		# fen / forested wetland
		lc_mod[ (lc_mod >= 41) & (lc_mod <= 95) & (canopy > 20) & (seak_mask == 1) ] = 11
		lc_mod[ (lc_mod >= 41) & (lc_mod <= 95) & (canopy <= 20) & (seak_mask == 1) ] = 12

		# harvested areas to upland
		logged = rasterio.open( input_paths[ 'logged' ] ).read_band( 1 )
		lc_mod[ (logged > 0) & (seak_mask == 1) ] = 10

		lc.write_band( 1, lc_mod )


	# mask the non-scak region <--> ugly...
	with rasterio.open( lc.name, mode='r+' ) as lc:
		lc_mod = lc.read_band( 1 )
		scak_mask = rasterio.open( input_paths[ 'scak_mask' ] ).read_band( 1 )
		# scak_mask = np.ma.masked_where( scak_mask==0, scak_mask )
		# cur_mask = np.ma.mask_or( scak_mask.mask, lc_mod_mask )
		# lc_mod = np.ma.masked_array( lc_mod.data, cur_mask )
		# lc.write_band( 1, lc_mod )

		# collapse initial undesired classes to noveg
		lc_mod = lc.read_band( 1 )
		lc_mod[ (lc_mod >= 0) & (lc_mod <= 31 ) & (scak_mask == 1) ] = 1 # potentially 0 later on

		# alder
		lc_mod[ (lc_mod == 41) & (scak_mask == 1) ] = 1		

		# white spruce
		lc_mod[ (lc_mod == 42) & (scak_mask == 1) ] = 2

		# deciduous
		lc_mod[ (lc_mod == 43) & (scak_mask == 1) ] = 3

		# shrub tundra
		lc_mod[ (lc_mod == 51) & (lc_mod == 52) & (lc_mod == 90) & (scak_mask == 1) ] = 4

		# heath
		lc_mod[ (lc_mod == 71) & (lc_mod == 72) & (scak_mask == 1) ] = 9

		# wetland tundra
		lc_mod[ (lc_mod == 95) & (scak_mask == 1) ] = 6

		lc.write_band( 1, lc_mod )
