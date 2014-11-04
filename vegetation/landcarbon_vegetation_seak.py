# # # # 
# reclassify NLCD 2001 Landcover for the LandCarbon project 
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

if __name__ == '__main__':
	import os, rasterio, fiona
	from rasterio import features
	from rasterio.warp import reproject, RESAMPLING
	import numpy as np
	
	# some initial setup
	version_num = 'v0_3'
	input_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/seak'
	output_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Output_Data'
	os.chdir( output_dir )
	
	input_paths = { 
			'lc01':os.path.join( input_dir, 'NLCD_land_cover_AKNPLCC.tif' ),
			'cp01':os.path.join( input_dir, 'NLCD_pct_canopy_AKNPLCC.tif' ),
			'logged':os.path.join( input_dir, 'AKNPLCC_2ndGrowth.tif' ),
			'seak_mask':os.path.join( input_dir, 'seak_aoi.tif' )			
	}

	# collapse initial undesired classes to noveg
	lc = rasterio.open( input_paths[ 'lc01' ] )
	lc_mod = lc.read_band( 1 )
	lc_mod[ (lc_mod >= 0) & (lc_mod <= 31 ) ] = 1 # potentially 0 later on

	# upland forest / fen
	canopy = rasterio.open( input_paths[ 'cp01' ] ).read_band( 1 )
	lc_mod[ (lc_mod == 42) & (canopy > 20) ] = 10
	lc_mod[ (lc_mod == 42) & (canopy <= 20) ] = 12

	# alder / shrubland
	lc_mod[ (lc_mod == 81) | (lc_mod == 82) ] = 15

	# fen / forested wetland
	lc_mod[ (lc_mod >= 41) & (lc_mod <= 95) & (canopy > 20) ] = 11
	lc_mod[ (lc_mod >= 41) & (lc_mod <= 95) & (canopy <= 20) ] = 12

	# harvested areas to upland
	logged = rasterio.open( input_paths[ 'logged' ] ).read_band( 1 )
	lc_mod[ (logged > 0) ] = 10

	meta = lc.meta
	meta.update( crs={'init': u'epsg:3338'} )

	output_filename = os.path.join( output_dir, 'landcarbon_model_vegetation_input_seak_2001_' + version_num + '.tif' )

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write_band( 1, lc_mod )
