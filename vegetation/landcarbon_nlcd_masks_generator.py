import rasterio, os, shutil, glob
import numpy as np

###
# REQUIRES GDAL 1.8.0 + 
#
# This is more of a roadmap to how we got to specific mask results than it is a script.
# Keep this in mind when reading the processes used to get to the needed masks.
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


## Now lets take the domain rasters for SEAK, SCAK and remove Canada 
ak_only_combined = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/combined_mask_alaska_only.tif' )
seak = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/seak_aoi.tif' )
scak = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/scak_aoi.tif' )
kodiak = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/kodiak_aoi.tif' )
# some file meta for the output rasters
meta = ak_only_combined.meta
meta.update( compress='lzw', nodata=None )

ak_arr = ak_only_combined.read_band( 1 )

# seak 
seak_arr = seak.read_band( 1 )
out_arr = np.zeros_like( seak_arr.data )

out_arr[ (ak_arr == 1) & (seak_arr == 1) ] = 1

output_filename = seak.name.replace( '.tif', '_akonly.tif' )
with rasterio.open( output_filename, 'w', **meta  ) as out:
	out.write_band( 1, out_arr )

del seak_arr, seak

# scak 
scak_arr = scak.read_band( 1 )
out_arr = np.zeros_like( scak_arr.data )

out_arr[ (ak_arr == 1) & (scak_arr == 1) ] = 1

output_filename = scak.name.replace( '.tif', '_akonly.tif' )
with rasterio.open( output_filename, 'w', **meta  ) as out:
	out.write_band( 1, out_arr )

del scak_arr, scak

# kodiak
kodiak_arr = kodiak.read_band( 1 ).filled()
output_filename = kodiak.name.replace( '.tif', '_akonly.tif' )
with rasterio.open( output_filename, 'w', **meta  ) as out:
    out.write_band( 1, kodiak_arr.filled() )

del kodiak_arr, kodiak


## Now lets take the domain rasters for SEAK, SCAK and remove Alaska
ak_only_combined = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/combined_mask_alaska_only.tif' )
seak = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/seak_aoi.tif' )
scak = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/scak_aoi.tif' )
kodiak = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/kodiak_aoi.tif' )
# some file meta for the output rasters
meta = ak_only_combined.meta
meta.update( compress='lzw', nodata=None )

ak_arr = ak_only_combined.read_band( 1 )

# seak 
seak_arr = seak.read_band( 1 )
out_arr = np.zeros_like( seak_arr.data )

out_arr[ (ak_arr == 0) & (seak_arr == 1) ] = 1

output_filename = seak.name.replace( '.tif', '_canadaonly.tif' )
with rasterio.open( output_filename, 'w', **meta  ) as out:
	out.write_band( 1, out_arr )

del seak_arr, seak

# scak 
scak_arr = scak.read_band( 1 )
out_arr = np.zeros_like( scak_arr.data )

out_arr[ (ak_arr == 0) & (scak_arr == 1) ] = 1

output_filename = scak.name.replace( '.tif', '_canadaonly.tif' )
with rasterio.open( output_filename, 'w', **meta  ) as out:
	out.write_band( 1, out_arr )

del scak_arr, scak


# # make a mask of the saltwater domain from NLCD -- class 11
nlcd = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/nlcd_2001_land_cover_maritime.tif' )
meta = nlcd.meta
meta.update( compress='lzw', nodata=None )
arr = nlcd.read_band( 1 ).data
arr[ arr != 11 ] = 0
arr[ arr == 11 ] = 1

output_filename = nlcd.name.replace( '.tif', '_saltwater_mask.tif' )
with rasterio.open( output_filename, 'w', **meta ) as out:
	out.write_band( 1, arr )


# make a 1km version in the domain of the large AKCanada Extent
output_resampled = out.name.replace( '.tif', '_akcan_1km.tif' )
if os.path.exists( output_resampled ):
	[ os.remove( i ) for i in glob.glob( output_resampled[:-3] + '*' ) ]

shutil.copyfile( template, output_resampled )
command = 'gdalwarp -r mode -multi -srcnodata None ' + out.name + ' ' + output_resampled
os.system( command )
