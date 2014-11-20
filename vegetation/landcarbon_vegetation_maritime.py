# # # # 
# reclassify NLCD 2001 Landcover for the LandCarbon project 
#	SEAK / SCAK / Kodiak Island Domains
# Author: Michael Lindgren (malindgren@alaska.edu) 
# # # # 

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
	hex = hex.lstrip( '#' )
	hlen = len( hex )
	rgb = [ int( hex[ i : i + hlen/3 ], 16 ) for i in range( 0, hlen, hlen/3 ) ]
	rgb.insert( len( rgb ) + 1, 1 )
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
	return { int( i.get( 'value' ) ) : tuple( hex_to_rgb( i.get( 'color' ) ) ) for i in tree.iter( tag='item' ) }

if __name__ == '__main__':
	import os, rasterio, fiona, shutil
	from rasterio import features
	from rasterio.warp import reproject, RESAMPLING
	import numpy as np

	# # some initial setup
	version_num = 'v0_4'
	input_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime'
	output_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Output_Data'
	qml_style = os.path.join( output_dir, 'qgis_styles','landcarbon_modeled_vegetation_2001_style.qml' )
	output_filename = os.path.join( output_dir, 'landcarbon_vegetation_modelinput_maritime_2001_' + version_num + '.tif' )

	os.chdir( output_dir )
	meta_updater = dict( driver='GTiff', dtype=rasterio.uint8, compress='lzw', crs={'init':'epsg:3338'}, count=1, nodata=255 )
	
	input_paths = { 
			'lc01':os.path.join( input_dir, 'nlcd_2001_land_cover_maritime.tif' ),
			'cp01':os.path.join( input_dir, 'nlcd_2001_canopy_cover_maritime.tif' ),
			'logged':os.path.join( input_dir, 'AKNPLCC_2ndGrowth.tif' ),
			'seak_mask':os.path.join( input_dir, 'seak_aoi_akonly.tif' ),
			'scak_mask':os.path.join( input_dir, 'scak_aoi_akonly.tif' ),
			'kodiak_mask':os.path.join( input_dir, 'kodiak_aoi_akonly.tif' )
	}

	# # open mask arrays
	seak_mask = rasterio.open( input_paths[ 'seak_mask' ] ).read_band( 1 )
	scak_mask = rasterio.open( input_paths[ 'scak_mask' ] ).read_band( 1 )
	kodiak_mask = rasterio.open( input_paths[ 'kodiak_mask' ] ).read_band( 1 )

	with rasterio.open( input_paths[ 'lc01' ], mode='r' ) as lc:
		lc_arr = lc.read_band( 1 )
		lc_arr.fill_value = 255
		lc_arr = lc_arr.filled()

		## ##  --- seak reclass --- ## ##
		# collapse initial undesired classes to noveg
		lc_arr = lc.read_band( 1 )
		lc_arr[ (lc_arr >= 0) & (lc_arr <= 31 ) ] = 1 # potentially 0 later on

		# upland forest / fen
		canopy = rasterio.open( input_paths[ 'cp01' ] ).read_band( 1 )
		lc_arr[ (lc_arr == 42) & (canopy > 20) & (seak_mask == 1) ] = 8
		lc_arr[ (lc_arr == 42) & (canopy <= 20) & (seak_mask == 1) ] = 10
		# alder / shrubland
		lc_arr[ (lc_arr == 81) | (lc_arr == 82) & (seak_mask == 1) ] = 11

		# fen / forested wetland
		lc_arr[ (lc_arr >= 41) & (lc_arr <= 95) & (canopy > 20) & (seak_mask == 1) ] = 9
		lc_arr[ (lc_arr >= 41) & (lc_arr <= 95) & (canopy <= 20) & (seak_mask == 1) ] = 10
		del canopy
		# harvested areas to upland
		logged = rasterio.open( input_paths[ 'logged' ] ).read_band( 1 )
		logged = logged.filled()
		lc_arr[ (logged == 1) & (seak_mask == 1) ] = 8
		del logged

		## ##  --- scak reclass --- ## ##
		# alder
		lc_arr[ (lc_arr == 41) & (scak_mask == 1) ] = 11	
		# white spruce
		lc_arr[ (lc_arr == 42) & (scak_mask == 1) ] = 2

		# deciduous
		lc_arr[ (lc_arr == 43) & (scak_mask == 1) ] = 3
		# shrub tundra
		lc_arr[ ( (lc_arr == 51) | (lc_arr == 52) | (lc_arr == 90) ) & (scak_mask == 1) ] = 4

		# heath
		lc_arr[ ( (lc_arr == 71) | (lc_arr == 72) ) & (scak_mask == 1) ] = 7
		# wetland tundra
		lc_arr[ (lc_arr == 95) & (scak_mask == 1) ] = 6

		# leftover 82 (Cultivated Crops) in the scak
		lc_arr[ (lc_arr == 82) & (scak_mask == 1) ] = 1

		## ##  --- kodiak reclass --- ## ##
		# deciduous
		lc_arr[ (lc_arr == 41) & (kodiak_mask == 1) ] = 3
		# upland forest
		lc_arr[ (lc_arr == 42) & (kodiak_mask == 1) ] = 8

		# shrub tundra
		lc_arr[ ( (lc_arr == 51) | (lc_arr == 90) ) & (kodiak_mask == 1) ] = 4
		# alder
		lc_arr[ (lc_arr == 52) & (kodiak_mask == 1) ] = 11

		# gramminoid tundra
		lc_arr[ (lc_arr == 71) & (kodiak_mask == 1) ] = 5
		# wetland tundra
		lc_arr[ ( (lc_arr == 72) | (lc_arr == 95) ) & (kodiak_mask == 1) ] = 6

		# set to not modeled any pixels that are outside the aoi masks
		combined_mask = np.sum( np.dstack( [seak_mask, scak_mask, kodiak_mask] ), axis=2 )
		del seak_mask, scak_mask, kodiak_mask
		combined_mask = combined_mask.astype( np.bool ) == False

		# fill nodata with 255 and set it as a masked array before writing to disk
		lc_arr = lc_arr.filled()
		lc_arr = np.ma.masked_array( lc_arr, combined_mask, fill_value=255, copy=True )

		# # write to disk
		meta = lc.meta
		meta.update( meta_updater )
		with rasterio.open( output_filename, mode='w', **meta ) as output:
			ctable = qml_to_ctable( qml_style ) # rasterio colormap dict r,g,b,a
			output.write_band( 1, lc_arr.filled() )
			output.write_colormap( 1, ctable )

