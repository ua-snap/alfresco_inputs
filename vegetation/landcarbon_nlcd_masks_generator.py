import rasterio, os, shutil
import numpy as np

###
# REQUIRES GDAL 1.8.0 + 
###

input_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Output_Data'
akcan = os.path.join( input_dir, 'alfresco_model_vegetation_input_2005.tif' )
template = os.path.join( input_dir, 'template_akcan_extent_1km.tif' )

# make a template alaska_canada empty raster for warping / resampling to
with rasterio.open( akcan ) as akcan_rst:
	meta = akcan_rst.meta
	meta.update( compress='lzw', nodata=None, dtype=np.uint8 )
	arr = akcan_rst.read_band( 1 ).data
	with rasterio.open( template, 'w', **meta ) as template_rst:
		template_rst.write_band( 1, np.zeros_like( arr ) )

# cleanup
del arr, akcan_rst, template_rst


## Maritime Domain Alaska Only 
combined_mask = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/combined_mask.tif' )
nlcd_mask = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/nlcd_2001_land_cover_maritime_mask.tif' )

combined_arr = combined_mask.read_band(1)
nlcd_arr = nlcd_mask.read_band(1)

nd = nlcd_arr.data
cd = combined_arr.data

new = np.zeros_like( cd )
new[ (cd == 1) & (nd == 1 ) ] = 1

new_mask = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/combined_mask_alaska_only.tif'
output_resampled = new_mask.replace( '.tif', '_akcan_1km.tif' )
meta = nlcd_mask.meta
meta.update( compress='lzw', nodata=None )

with rasterio.open( new_mask, 'w', **meta ) as out:
	out.write_band( 1, new )

# make a 1km version in the domain of the large AKCanada Extent
if os.path.exists( output_resampled ):
	[ os.remove( i ) for i in glob.glob( output_resampled[:-3] + '*' ) ]

shutil.copyfile( template, output_resampled )
command = 'gdalwarp -r mode -multi -srcnodata None ' + new_mask + ' ' + output_resampled
os.system( command )


## Maritime Domain Canada Only
new = np.zeros_like( cd )
new[ (cd == 1) & (nd == 0 ) ] = 1

new_mask = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/combined_mask_canada_only.tif'
output_resampled = new_mask.replace( '.tif', '_akcan_1km.tif' )
meta = nlcd_mask.meta
meta.update( compress='lzw', nodata=None )
with rasterio.open( new_mask, 'w', **meta ) as out:
	out.write_band( 1, new )

# make a 1km version in the domain of the large AKCanada Extent
if os.path.exists( output_resampled ):
	[ os.remove( i ) for i in glob.glob( output_resampled[:-3] + '*' ) ]

shutil.copyfile( template, output_resampled )
command = 'gdalwarp -r mode -multi -srcnodata None ' + new_mask + ' ' + output_resampled
os.system( command )
