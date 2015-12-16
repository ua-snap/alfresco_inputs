library(sp)
library(maptools)
library(rgdal)
library(raster)


# read in the stack of the 12 months of the variable to downscale to
# read in the DAT file from NEW et al 2002
t <- read.table("/Data/Base_Data/Climate/World/CRU_grids/CRU_TS20/grid_10min_tmp.dat")
# set the column names to the lon/lat and months
colnames(t) <- c("lat","lon","jan","feb","mar","apr","jun","jul","aug","sep","oct","nov","dec")
cols <- c("jan","feb","mar","apr","jun","jul","aug","sep","oct","nov","dec")

template <- raster("/big_storage/malindgren/AIEM/RHUM/tas_cru_10min_1961_1990/grid_latlong_template/hur_cru_10min_01_1961_1990.tif")

#snap.extent <- readOGR("/workspace/UA/malindgren/C_Michael_Data/DATA/Shapefiles/AKCanada_Downscaled_extent/AK_Canada_DissolvedExtent_wNWT.shp", "AK_Canada_DissolvedExtent_wNWT")
count=0
# now lets loop through the data columns and turn them into rasters
for(i in 3:14){
	count=count+1

	if(nchar(count) < 2){ count2 = paste("0",count, sep="") } else { count2=count }

	# turn the lon lat into a points object
	t.points <- SpatialPoints(coords=t[,2:1])
	
	# then turn that into spatialpoints dataframe
	t.spdf <- SpatialPointsDataFrame(coords=t.points, data=data.frame(t[,i]))

	names(t.spdf) <- "tmp"
	
	# clip that point layer with a polygon of the extent that fits SNAPs full extent (in wgs84)
	# writeOGR(t.spdf, dsn="/big_storage/malindgren/AIEM/RHUM/tas_cru_10min_1961_1990/point_shape/", layer=paste("cru_10min_tmp_1961_1990_",count, sep=""), driver="ESRI Shapefile")

	# Now turn that spdf into a raster object
	r <- rasterize(t.spdf, template, field="tmp")

	writeRaster(r, filename=paste("/big_storage/malindgren/AIEM/RHUM/tas_cru_10min_1961_1990/point_shape/", "cru_10min_tmp_1961_1990_",count2, ".tif",sep=""))

}
