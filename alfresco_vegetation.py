# # # # 
# reclassify input NALCMS LandCover Map to classes useful in the ALRESCO Fire Model
#  using pre-processed ancillary layers.
# # # #

# NOTES:  
# a little neighborhood math for the queens case
# ul = (-1,-1)
# uc = (-1, 0)
# ur = (-1,+1)
# l = (0,-1)
# r = (0,+1)
# ll = (+1,-1)
# lc = (+1,0)
# lr = (+1,+1)


def replace_erroneous_treeline( lc2d, tl2d ):
	''' replace values based on neighbors of offenders and a where condition 
		arguments:
			lc2d = 2d landcover array with erroneous values within
			tl2d = 2d treeline boolean array to subset the landcover array
		returns:
			2d numpy array with the erroneous values removed

	'''
	ind = np.where( (((lc2d == 1) | (lc2d == 2)) & (tl2d == 1)) ) # this is hardwired...
	ind_zip = zip( *ind )
	print len( ind_zip )
	if len( ind_zip ) == 0:
		return lc2d
	else:
		index_groups = [ [(i-1,j-1), 
						  (i-1, j+0), 
						  (i-1,j+1), 
						  (i+0,j-1), 
						  (i+0,j+1), 
						  (i+1,j-1), 
						  (i+1,j+0), 
						  (i+1,j+1)]
						for i,j in ind_zip ]
		new_ind = []
		for count, group in enumerate( index_groups ):
			cols = np.array( [j for i,j in group] )
			rows = np.array( [i for i,j in group] )
			vals = lc2d[ ( rows, cols ) ]
			vals = vals[ (vals!=1) & (vals!=2) ] # (vals > 0) & 
			if len(vals) != 0:
				uniques, counts = np.unique( vals, return_counts = True )
				new_val = uniques[ np.argmax( counts ) ]
				if new_val not in [1,2]:
					lc2d[ ind_zip[ count ] ] = new_val
		return replace_erroneous_treeline( lc2d, tl2d )


if __name__ == '__main__':
	import rasterio, fiona, os, sys
	import numpy as np
	import scipy as sp

	curdir = os.getcwd()
	os.chdir( '/workspace/Shared/Tech_Projects/AK_LandCarbon/project_data/CODE' )
	from geolib_snap import reclassify

	# this should be something passed in at the command line
	gs_value = 6.5
	output_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Output_Data'
	output_veg = os.path.join( output_dir, 'alfresco_model_vegetation_input_2005.tif' )
	input_paths = {
		'lc05':"/workspace/UA/malindgren/projects/NALCMS_Veg_reClass/August2012_FINALversion/ALFRESCO_VegMap_Ancillary/na_landcover_2005_1km_MASTER.tif",
		'north_south':"/workspace/UA/malindgren/projects/NALCMS_Veg_reClass/August2012_FINALversion/ALFRESCO_VegMap_Ancillary/AKCanada_1km_NorthSouth_FlatWater_999_MASTER.tif",
		'mask':"/workspace/UA/malindgren/projects/NALCMS_Veg_reClass/August2012_FINALversion/ALFRESCO_VegMap_Ancillary/mask_for_finalization_alfresco_VegMap.tif",
		'gs_temp':"/workspace/UA/malindgren/projects/NALCMS_Veg_reClass/August2012_FINALversion/ALFRESCO_VegMap_Ancillary/AKCanada_gs_temp_mean_MJJAS_1961_1990_climatology_1km_bilinear_MASTER.tif",
		'coast_spruce_bog':"/workspace/UA/malindgren/projects/NALCMS_Veg_reClass/August2012_FINALversion/ALFRESCO_VegMap_Ancillary/Coastal_vs_Woody_wetlands_MASTER.tif",
		'treeline':"/workspace/UA/malindgren/projects/NALCMS_Veg_reClass/August2012_FINALversion/ALFRESCO_VegMap_Ancillary/CAVM_treeline_AKCanada_1km_commonExtent_MASTER.tif",
		'NoPac':"/workspace/UA/malindgren/projects/NALCMS_Veg_reClass/August2012_FINALversion/ALFRESCO_VegMap_Ancillary/ALFRESCO_NorthPacMaritime_forVegMap.tif"
	}

	# collapse classes in the input landcover raster
	lc = rasterio.open( input_paths[ 'lc05' ] )
	reclass_list = [[15,20,0],[128,129,0],[1,3,9],[5,7,3],[11,12,4]]
	output_filename = os.path.join( output_dir, 'alfresco_vegetation_step1.tif' )
	lc_mod = reclassify( lc, reclass_list, output_filename, band=1, creation_options=dict() )
	lc_mod = lc_mod.read_band( 1 )
	lc_mod.fill_value = 9999

	# convert wetland to spruce bog and wetland
	coast_spruce_bog = rasterio.open( input_paths[ 'coast_spruce_bog' ] ).read_band( 1 )
	coast_spruce_bog.fill_value = 9999
	lc_mod[ (lc_mod == 14) & (coast_spruce_bog == 2) ] = 9
	lc_mod[ (lc_mod == 14) & (coast_spruce_bog != 2) ] = 20

	# coastal wetland class to WETLAND TUNDRA or NO VEG based on the gs_temp (Average Growing Season Temperature) values
	treeline = rasterio.open( input_paths[ 'treeline' ] ).read_band(1).astype( np.int32 )
	gs_temp = rasterio.open( input_paths[ 'gs_temp' ] ).read_band(1)
	treeline.fill_value = 9999
	gs_temp.fill_value = 9999.0
	lc_mod[ (lc_mod == 20) & (gs_temp < gs_value) & (treeline == 1) ] = 9
	lc_mod[ (lc_mod == 20) ] = 0

	# turn the placeholder class 8 (Temperate or sub-polar shrubland) into DECIDUOUS or SHRUB TUNDRA
	# NOTE: may be best to do the treeline query here for these weird values <- (treeline == 1) etc
	lc_mod[ (lc_mod == 8) & (gs_temp < gs_value) ] = 4
	lc_mod[ (lc_mod == 8) & (gs_temp >= gs_value) ] = 3

	# Reclass Sub-polar or polar grassland-lichen-moss as GRAMMINOID TUNDRA
	lc_mod[ (lc_mod == 10) | (lc_mod == 12) ] = 5

	# Take the SPRUCE placeholder class and parse it out in to WHITE / BLACK.
	# White Spruce = SPRUCE class & on a South-ish Facing slope
	north_south = rasterio.open( input_paths[ 'north_south' ] ).read_band(1)
	north_south.fill_value = 9999
	lc_mod[ (lc_mod == 9) & (north_south == 1) ] = 2
	lc_mod[ (lc_mod == 9) & (north_south == 2) ] = 1
	lc_mod[ (lc_mod == 9) ] = 1 # convert low lying spruce leftover to black spruce

	# reclassify erroneous spruce pixels with the most common values of its neighbors
	lc_mod = replace_erroneous_treeline( lc_mod, treeline )
	
	# temperate rainforest (seak) region delineation
	temperate_rainforest = rasterio.open( input_paths[ 'NoPac' ] ).read_band( 1 )
	lc_mod[ (lc_mod > 0) & (temperate_rainforest == 1)  ] = 8

	# conversion of the barren-lichen moss
	lc_mod[ lc_mod==13 ] = 7

	# convert oob to 255
	mask = rasterio.open( input_paths[ 'mask' ] ).read_band( 1 )
	lc_mod[ mask == 1 ] = 255

	# change the type (byte) and write it out
	meta = lc.meta.copy()
	meta.update( dtype=rasterio.uint8 )

	with rasterio.open( output_veg, 'w', **meta ) as out:
		out.write_band( 1, lc_mod.astype( rasterio.uint8 ))

