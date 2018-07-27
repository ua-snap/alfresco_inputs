# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# rasterize the Alaska Fire Perimeters from the AICC derived GDB and the Canada Fire Perimeters from the NFDB
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

import rasterio, os, glob
from rasterio.features import rasterize
import geopandas as gpd
import numpy as np

# # input variables
firehistory_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Fire/FireHistory/raw_downloaded/FireHistoryPerimeters_1940_2017.gdb'
layername = 'FireHistoryPerimeters'
out_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Fire/FireHistory/AK_CAN_1km'
template_fn = '/Data/Base_Data/ALFRESCO/ALFRESCO_Master_Dataset_v2_1/ALFRESCO_Model_Input_Datasets/AK_CAN_Inputs/FirePrep/template/alf_v1_mask_1km.tif'
firehistory_fn_canada = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Fire/FireHistory/raw_downloaded/NFDB_poly_20171106.shp'

# # # # # ALASKA DATA # # # # 
df = gpd.read_file( firehistory_fn, layer=layername )
df = df.to_crs( epsg=3338 )

# # make a template of the meta
with rasterio.open( template_fn ) as rst:
	meta = rst.meta.copy()
	meta.update( compress='lzw', nodata=-9999 )
	shape = rst.shape
	mask = (rst.read_masks(1) > 0).astype(int)

# # # # # CANADA DATA # # # # 
can_df = gpd.read_file( firehistory_fn_canada )
can_df = can_df.to_crs( epsg=3338 )
can_df = can_df[ can_df.YEAR != -9999 ]

can = pd.DataFrame({'geometry': can_df.geometry, 'year':can_df.YEAR.astype(int) })
ak = pd.DataFrame({'geometry': df.geometry, 'year':df.FIREYEAR.astype(int) })
akcan = gpd.GeoDataFrame( pd.concat([ak, can]), geometry='geometry', crs={'init':'epsg:3338'} )

for group, sub_df in akcan.groupby( 'year' ):
	print( group )
	geoms = sub_df.geometry.copy()
	out_rst = rasterize( geoms, out_shape=shape, fill=0, transform=meta['transform'], all_touched=False, default_value=1, dtype=np.float32 )
	out_rst[ mask == 0 ] = -9999

	output_filename = os.path.join( out_path, 'ALF_AK_CAN_FireHistory_{}.tif'.format( year ) )
	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( out_rst, 1 )
