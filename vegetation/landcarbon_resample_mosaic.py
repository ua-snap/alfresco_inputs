with rasterio.open( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Vegetation/Input_Data/maritime/land_sea_nlcd_2.tif', 'w', **meta ) as out:
	out.write_band( 1, arr.astype( np.uint8 ) )




command = 'gdalwarp -tr 1000 1000 -r mode -srcnodata None -dstnodata None -multi -co "COMPRESS=LZW" '+ output_lc.name + ' ' + output_filename
os.system( command )
