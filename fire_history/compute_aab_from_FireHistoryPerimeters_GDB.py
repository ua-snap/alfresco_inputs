# compute annual area burned in 
import rasterio, os, glob
from rasterio.features import rasterize
import pandas as pd
import geopandas as gpd
import numpy as np

firehistory_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Fire/FireHistory/raw_downloaded/FireHistoryPerimeters_1940_2017.gdb'
layername = 'FireHistoryPerimeters'
output_filename = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Fire/FireHistory/Alaska_AnnualAreaBurned_km2_1940-2017.csv'

# read the data and make sure it is in 3338 --> units: meters
df = gpd.read_file( firehistory_fn, layer=layername )
df = df.to_crs( epsg=3338 )

# groupby year and get the total area burned for each / convert to km2 / round it / convert to int / output as dataframe
aab = df.geometry.groupby( df['FIREYEAR'] ).apply( lambda x: x.apply( lambda y: y.area).sum() ).div( 1000000 ).round( 0 ).astype( np.int32 ).to_frame( name='AnnualAreaBurned' )
aab.to_csv( output_filename, sep=',' ) # write file to csv

# to read the csv file back in with pandas -- then you can plot and do all kinds of stuff with this DataFrame object
df = pd.read_csv( output_filename, index_col=0 )



# # # # # CANADA DATA -- Still stuck in 2016...

firehistory_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Fire/FireHistory/raw_downloaded/NFDB_poly_20171106.shp'
can_df = gpd.read_file( firehistory_fn )
can_df = can_df.to_crs( epsg=3338 )
aab_canada = df.geometry.groupby( df['YEAR'] ).apply( lambda x: x.apply( lambda y: y.area).sum() ).div( 1000000 ).round( 0 ).astype( np.int32 ).to_frame( name='AnnualAreaBurned' )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
