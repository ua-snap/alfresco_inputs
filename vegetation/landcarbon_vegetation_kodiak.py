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

	# import local library of functions
	os.chdir( '/workspace/Shared/Tech_Projects/AK_LandCarbon/project_data/CODE' )
	from geolib_snap import reclassify

	# some initial setup
	version_num = 'v0_4'
	input_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/kodiak'
	output_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Output_Data'
	os.chdir( output_dir )
	
	# # THIS NEEDS UPDATING TO THE KODIAK MAP NEEDS!
	input_paths = { 
			'lc01':os.path.join( input_dir, 'ak_nlcd_2001_land_cover_3130_KodiakIsland_3338.tif' ),
			'kodiak_mask':os.path.join( input_dir, 'kodiak_aoi_box.tif' )
	}

	# reclassify
	lc = rasterio.open( input_paths[ 'lc01' ] )

	output_filename = os.path.join( output_dir, 'landcarbon_model_vegetation_input_kodiak_2001_' + version_num + '.tif' )
	reclass_list = [[11,32,1],[41,42,3],[42,43,10],[51,52,4],[52,53,15],[71,72,5],[72,73,6],[90,91,4],[95,96,6]]
	out_rst = reclassify( lc, reclass_list, output_filename, band=1, creation_options={'compress':'lzw'} )
	out_rst.close()
