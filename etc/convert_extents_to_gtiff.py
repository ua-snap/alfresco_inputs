# # # # # #
# convert the shapefile extents data to boolean geotiffs
# # # # # #

import glob, os, fiona, rasterio
from rasterio import features


l = glob.glob( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/extents/*.shp' )

shapes = [ [ (j['geometry'], 1) for j in fiona.open(i)] for i in l ]

final = [ features.rasterize(i, out_shape=lc.shape, transform=lc.transform, fill=0 ) for i in shapes ]

for i,j in zip(final, l):
    rst = rasterio.open( j.replace('.shp','.tif'), 'w', **meta )
    rst.write_band(1, i)
    rst.close()
    