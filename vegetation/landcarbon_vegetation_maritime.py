# # # # 
# reclassify NLCD 2001 Landcover for the LandCarbon project 
#	SEAK / SCAK / Kodiak Island Domains
# # # # 

# Dev Notes:
# ---
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


# # # OUTPUT LEGEND # # # #
#
#	1 - not modeled
#	2 - white spruce
#	3 - deciduous
#	4 - shrub tundra
#	5 - gramminoid tundra
#	6 - wetland tundra
#	7 - heath
#	8 - upland forest
#	9 - forested wetland
#	10 - fen
#	11 - alder
#	255 - out of bounds
# # # # # # # # # # # # # # # 


if __name__ == '__main__':
	import os, rasterio, fiona, shutil
	from rasterio import features
	from rasterio.warp import reproject, RESAMPLING
	import numpy as np

	# # some initial setup
	version_num = 'v0_4'
	input_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime'
	output_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Output_Data'
	output_filename = os.path.join( output_dir, 'landcarbon_vegetation_modelinput_maritime_2001_' + version_num + '.tif' )

	os.chdir( output_dir )
	meta_updater = dict( driver='GTiff', dtype=rasterio.uint8, compress='lzw', crs={'init':'epsg:3338'}, count=1, nodata=None )
	
	input_paths = { 
			'lc01':os.path.join( input_dir, 'nlcd_2001_land_cover_maritime.tif' ),
			'cp01':os.path.join( input_dir, 'nlcd_2001_canopy_cover_maritime.tif' ),
			'logged':os.path.join( input_dir, 'AKNPLCC_2ndGrowth.tif' ),
			'seak_mask':os.path.join( input_dir, 'seak_aoi.tif' ),
			'scak_mask':os.path.join( input_dir, 'scak_aoi.tif' ),
			'kodiak_mask':os.path.join( input_dir, 'kodiak_aoi.tif' )
	}

	# # open mask arrays
	scak_mask = rasterio.open( input_paths[ 'scak_mask' ] ).read_band( 1 )
	seak_mask = rasterio.open( input_paths[ 'seak_mask' ] ).read_band( 1 )
	kodiak_mask = rasterio.open( input_paths[ 'kodiak_mask' ] ).read_band( 1 )

	with rasterio.open( input_paths[ 'lc01' ], mode='r' ) as lc:
		lc_arr = lc.read_band( 1 )
		
		## ## seak reclass ## ##
		# collapse initial undesired classes to noveg
		lc_arr = lc.read_band( 1 )
		lc_arr[ (lc_arr >= 0) & (lc_arr <= 31 ) ] = 0 # potentially 0 later on

		# upland forest / fen
		canopy = rasterio.open( input_paths[ 'cp01' ] ).read_band( 1 )
		lc_arr[ (lc_arr == 42) & (canopy > 20) & (seak_mask == 1) ] = 8 #--#
		lc_arr[ (lc_arr == 42) & (canopy <= 20) & (seak_mask == 1) ] = 10 #--#
		# alder / shrubland
		lc_arr[ (lc_arr == 81) | (lc_arr == 82) & (seak_mask == 1) ] = 11 #--#

		# fen / forested wetland
		lc_arr[ (lc_arr >= 41) & (lc_arr <= 95) & (canopy > 20) & (seak_mask == 1) ] = 9 #--#
		lc_arr[ (lc_arr >= 41) & (lc_arr <= 95) & (canopy <= 20) & (seak_mask == 1) ] = 10 #--#
		# harvested areas to upland
		logged = rasterio.open( input_paths[ 'logged' ] ).read_band( 1 )
		lc_arr[ (logged == 1) & (seak_mask == 1) ] = 8 #--#

		## ## scak reclass ## ##
		# alder
		lc_arr[ (lc_arr == 41) & (scak_mask == 1) ] = 11 #--#	
		# white spruce
		lc_arr[ (lc_arr == 42) & (scak_mask == 1) ] = 2 #--#

		# deciduous
		lc_arr[ (lc_arr == 43) & (scak_mask == 1) ] = 3 #--#
		# shrub tundra
		lc_arr[ ( (lc_arr == 51) | (lc_arr == 52) | (lc_arr == 90) ) & (scak_mask == 1) ] = 4 #--#

		# heath
		lc_arr[ ( (lc_arr == 71) | (lc_arr == 72) ) & (scak_mask == 1) ] = 7 #--#
		# wetland tundra
		lc_arr[ (lc_arr == 95) & (scak_mask == 1) ] = 6 #--#

		## ## kodiak reclass ## ##
		# deciduous
		lc_arr[ (lc_arr == 41) & (kodiak_mask == 1) ] = 3 #--#
		# upland forest
		lc_arr[ (lc_arr == 42) & (kodiak_mask == 1) ] = 8 #--#

		# shrub tundra
		lc_arr[ ( (lc_arr == 51) | (lc_arr == 90) ) & (kodiak_mask == 1) ] = 4 #--#
		# alder
		lc_arr[ (lc_arr == 52) & (kodiak_mask == 1) ] = 11 #--#

		# gramminoid tundra
		lc_arr[ (lc_arr == 71) & (kodiak_mask == 1) ] = 5 #--#
		# wetland tundra
		lc_arr[ ( (lc_arr == 72) | (lc_arr == 95) ) & (kodiak_mask == 1) ] = 6 #--#

		# # write to disk
		meta = lc.meta
		meta.update( meta_updater )
		with rasterio.open( output_filename, mode='w', **meta ) as output:
			output.write_band( 1, lc_arr )
			# output.write_colormap( 1, ... ) # not implemented yet
