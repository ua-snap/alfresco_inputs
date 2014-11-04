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

if __name__ == '__main__':
	import os, rasterio, fiona
	from rasterio import features
	from rasterio.warp import reproject, RESAMPLING
	import numpy as np

	# import local library of functions
	os.chdir( '/workspace/Shared/Tech_Projects/AK_LandCarbon/project_data/CODE' )
	from geolib_snap import reclassify

	# some initial setup
	version_num = 'v0_3'
	input_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/seak'
	output_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Output_Data'
	os.chdir( output_dir )
	meta_updater = dict( driver='GTiff', dtype=rasterio.int16, compress='lzw', crs={'init':'EPSG:3338'}, count=1, nodata=None )
	
	input_paths = { 
			'lc01':os.path.join( input_dir, 'NLCD_land_cover_AKNPLCC.tif' ),
			'cp01':os.path.join( input_dir, 'NLCD_pct_canopy_AKNPLCC.tif' ),
			'logged':os.path.join( input_dir, 'AKNPLCC_2ndGrowth.tif' ),
			'scak_mask':os.path.join( input_dir, 'scak_aoi.tif' )
	}

	# reclassify
	lc = rasterio.open( input_paths[ 'lc01' ] )

	output_filename = os.path.join( output_dir, 'landcarbon_model_vegetation_input_scak_2001_' + version_num + '.tif' )
	reclass_list = [ [0,32,1],[41,42,15],[42,43,2],[43,44,1],[51,53,4],[71,73,9],[90,91,4],[95,96,6] ]
	output_rst = reclassify( lc, reclass_list, output_filename, band=1, creation_options={'compress'='lzw'} )
	output_rst.close()
