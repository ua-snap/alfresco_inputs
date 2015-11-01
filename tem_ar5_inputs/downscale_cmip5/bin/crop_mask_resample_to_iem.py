def resample_to_1km( x, template_raster_obj ):
	'''
	template_raster_obj should be a mask in in the res/extent/origin/crs of the 
	existing TEM IEM products.
	'''
	import rasterio, os
	from rasterio.warp import RESAMPLING, reproject
	import numpy as np

	fn = os.path.basename( x )
	fn_split = fn.split( '.' )[0].split( '_' ) 
	if '_cru_' in fn:
		output_path = os.path.dirname( x ).replace( '/cru_ts31/', '/IEM/cru_ts31/' ) # hardwired!
		fn_parts = ['variable', 'metric', 'model_1', 'model_2', 'kind', 'month', 'year']
		fn_dict = dict( zip( fn_parts, fn_split ) )
		fn_dict.update( scenario='historical', model='cru_ts31' )
	else:
		output_path = os.path.dirname( x ).replace( '/ar5/', '/IEM/ar5/' ) # hardwired!
		fn_parts = ['variable', 'metric', 'model', 'scenario', 'ensemble', 'month', 'year']
		fn_dict = dict( zip( fn_parts, fn_split ) )

	if not os.path.exists( output_path ):
		os.makedirs( output_path )

	fn_switch = { 'cld':'_'.join([ 'cld','mean','pct','iem',fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif',
		'vap':'_'.join(['vap','mean','hPa','iem', fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif', 
		'tas':'_'.join(['tas','mean','C','iem',fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif',
		'hur':'_'.join(['hur','mean','pct','iem',fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif' }

	output_filename = os.path.join( output_path, fn_switch[ fn_dict[ 'variable' ] ] )

	template_arr = template_raster_obj.read( 1 )
	template_meta = template_raster_obj.meta
	template_meta.update( compress='lzw' )
	
	rst = rasterio.open( x )
	# window = rst.window( *template_raster_obj.bounds )
	# window_transform = rst.window_transform( window )
	rst_arr = rst.read( 1 ) # window=window )
	
	if 'transform' in template_meta.keys():
		template_meta.pop( 'transform' )

	output_arr = np.empty_like( template_arr )

	mask = template_raster_obj.read_mask( 1 )
	mask[ mask != 0] = 1
	meta = template_raster_obj.meta
	meta.update( compress='lzw', dtype=np.uint8, nodata=0 )
	meta.pop( 'transform' )

	src_crs = {'init':'epsg:3338'}
	dst_crs = {'init':'epsg:3338'}
	reproject( rst_arr, output_arr, src_transform=rst.affine, src_crs=rst.crs, src_nodata=rst.nodata, \
			dst_transform=template_meta['affine'], dst_crs=template_meta['crs'],\
			dst_nodata=None, resampling=RESAMPLING.cubic_spline, num_threads=1 )

	with rasterio.open( output_filename, 'w', **template_meta ) as out:
		output_arr[ mask == 0 ] = template_raster_obj.nodata
		out.write( output_arr, 1 )
	return output_filename


if __name__ == '__main__':
	# some setup:
	template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/extents/vap_mean_hPa_iem_cccma_cgcm3_1_sresa1b_01_2001.tif'
	template_raster_obj = rasterio.open( template_raster_fn )
	resample_to_1km( x, template_raster_obj )


	with rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/extents/iem_mask_1km.tif', 'w', **meta ) as out:
		out.write( mask.astype(np.uint8), 1 )

	x = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_ts31_old/cld/downscaled/cld_pct_cru_ts31_downscaled_03_2001.tif'
	# x = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/ar5/IPSL-CM5A-LR/cld/downscaled/cld_metric_IPSL-CM5A-LR_rcp85_r1i1p1_12_2087.tif'



# def standardized_fn_to_vars( fn ):
# 	''' take a filename string following the convention for this downscaling and break into parts and return a dict'''
# 	fn = os.path.basename( fn )
# 	fn_split = fn.split( '.' )[0].split( '_' ) 
# 	if '_cru_' in fn:
# 		fn_parts = ['variable', 'metric', 'model_1', 'model_2', 'kind', 'month', 'year']
# 		fn_dict = dict( zip( fn_parts, fn_split ) )
# 		fn_dict.update( scenario='historical', model='cru_ts31' )
# 		# name_convention = [ 'variable', 'metric', 'model', 'scenario', 'experiment', 'begin_time', 'end_time' ]
# 	else:
# 		fn_parts = ['variable', 'metric', 'model', 'scenario', 'ensemble', 'month', 'year']
# 		fn_dict = dict( zip( fn_parts, fn_split ) )


# 	fn_switch = { 'cld':'_'.join([ 'cld','mean','pct','iem',fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif',
# 		'vap':'_'.join(['vap','mean','hPa','iem', fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif', 
# 		'tas':'_'.join(['tas','mean','C','iem',fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif',
# 		'hur':'_'.join(['hur','mean','pct','iem',fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif' }

# 	output_filename = os.path.join( output_path, fn_switch[ fn_dict[ 'variable' ] ] )

# 	fn_list = fn.split( '.' )[0].split( '_' )
# 	return { i:j for i,j in zip( name_convention, fn_list )}
