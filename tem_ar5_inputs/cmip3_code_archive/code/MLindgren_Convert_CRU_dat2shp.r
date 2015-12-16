library(sp)
library(maptools)
library(rgdal)


# read in the stack of the 12 months of the variable to downscale to
# read in the DAT file from NEW et al 2002
t <- read.table("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/extracted_files/grid_10min_reh.dat")
# set the column names to the lon/lat and months
colnames(t) <- c("lat","lon","jan","feb","mar","apr","jun","jul","aug","sep","oct","nov","dec")
cols <- c("jan","feb","mar","apr","jun","jul","aug","sep","oct","nov","dec")

#snap.extent <- readOGR("/workspace/UA/malindgren/C_Michael_Data/DATA/Shapefiles/AKCanada_Downscaled_extent/AK_Canada_DissolvedExtent_wNWT.shp", "AK_Canada_DissolvedExtent_wNWT")
count=0
# now lets loop through the data columns and turn them into rasters
for(i in 3:14){
	count=count+1
	# turn the lon lat into a points object
	t.points <- SpatialPoints(coords=t[,2:1])
	
	# then turn that into spatialpoints dataframe
	t.spdf <- SpatialPointsDataFrame(coords=t.points, data=data.frame(t[,i]))

	names(t.spdf) <- "rhum"
	
	# clip that point layer with a polygon of the extent that fits SNAPs full extent (in wgs84)
	writeOGR(t.spdf, dsn="/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/extracted_files/", layer=paste("cru_10min_rhum_1961-1990_",count, sep=""), driver="ESRI Shapefile")

	# Now turn that spdf into a raster object

}
