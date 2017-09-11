import os, shutil, rasterio

mask = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/alaska_canada/alfresco_vegetation_mask.tif'
# nalcms = '/workspace/Shared/Users/malindgren/SCOTT_August_2015/conifer_decid_ratios/nalcms10/NA_LandCover_2010_25haMMU.tif'
nalcms = '/workspace/Shared/Users/malindgren/SCOTT_August_2015/conifer_decid_ratios/nalcms10/NA_LandCover_2010_25haMMU_3338_arc.tif'

new_name = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/alaska_canada/NA_LandCover_2010_25haMMU_akcan_3.tif'
shutil.copy( mask, new_name )
# output_filename = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/alaska_canada/NA_LandCover_2010_25haMMU_test.tif'
arr = rasterio.open( mask ).read( 1 )
meta = rasterio.open( mask ).meta
meta.update( compress='lzw', crs={ 'init':'epsg:3338' } )
with rasterio.open( new_name, 'w', **meta ) as out:
	out.write( arr, 1 )

command = 'gdalwarp -multi ' + nalcms + ' ' + new_name
os.system( command )

rst = rasterio.open( new_name )
mask_arr = rasterio.open( mask ).read( 1 )
arr = rst.read( 1 )
arr[ mask_arr == 0 ] = 255
meta = rst.meta

meta.update( nodata=255 )
with rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/alaska_canada/NA_LandCover_2010_25haMMU_akcan_masked.tif', 'w', **meta ) as out:
	out.write( arr, 1 )

os.remove( new_name )




