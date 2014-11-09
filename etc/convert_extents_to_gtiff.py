# # # # # #
# convert the shapefile extents data to boolean geotiffs
# # # # # #

import glob, os, fiona, rasterio
from rasterio import features

lc = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/nlcd_2001_land_cover_maritime.tif' )
meta = lc.meta
meta.update( compress='lzw', crs={'init':'epsg:3338'} )

l = glob.glob( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/shapefile_extents/*.shp' )

shapes = [ [ (j['geometry'], 1) for j in fiona.open(i)] for i in l ]

final = [ features.rasterize(i, out_shape=lc.shape, transform=lc.transform, fill=0 ) for i in shapes ]

for i,j in zip(final, l):
	rst = rasterio.open( j.replace('.shp','.tif'), 'w', **meta )
	rst.write_band(1, i)
	rst.close()


# # # # extend the extent of the canopy raster

import glob, os, fiona, rasterio
from rasterio import features

lc = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/nlcd_2001_land_cover_maritime.tif' )

arr = lc.read_band( 1 )
dat = arr.data
dat[:] = 199

meta = lc.meta
meta.update( compress='lzw', crs={'init':'epsg:3338'} )

with rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/nlcd_maritime_full_bbox.tif', 'w', **meta ) as out:
	out.write_band( 1, dat )

os.system( 'gdal_merge.py -co COMPRESS=LZW -of GTiff -o /workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/nlcd_2001_canopy_cover_maritime.tif /workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/nlcd_maritime_full_bbox.tif /workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/NLCD_pct_canopy_AKNPLCC.tif' )

# # # STEP 2
with rasterio.open('/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/nlcd_2001_canopy_cover_maritime.tif', 'r+' ) as rst:
	arr = rst.read_band( 1 )
	arr[arr == 199 ] = 255
	rst.write_band( 1, arr )


# # # # extend the extent of the second growth raster by way of the shapefile input 
import glob, os, fiona, rasterio
from rasterio import features

lc = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/nlcd_2001_land_cover_maritime.tif' )
meta = lc.meta
meta.update( compress='lzw', crs={'init':'epsg:3338'} )
shp_name = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/AKNPLCC_2ndGrowth.shp'
shp = fiona.open( shp_name )
shapes = [ (j['geometry'], 1) for j in shp ]
final = features.rasterize(shapes, out_shape=lc.shape, transform=lc.transform, fill=0 )
with rasterio.open( shp_name.replace('.shp','.tif'), 'w', **meta ) as rst:
	rst.write_band( 1, final )

# # # MORE STUFF RELATED TO THE MASKING 

os.chdir( '/workspace/Shared/Tech_Projects/AK_LandCarbon/project_data/CODE' )
from geolib_snap.py import overlay_cover

full_aoi = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/nlcd_2001_land_cover_maritime_boolean.tif' ).read_band( 1 )
seak_aoi = rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/seak_aoi.tif' ).read_band( 1 )



overlay_cover( full_aoi, seak_aoi )


