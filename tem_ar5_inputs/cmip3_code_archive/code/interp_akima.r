d.begin <- date()

# read in the spatial packages needed in this downscaling
require(raster)
require(akima)
require(maptools)
require(rgdal)
require(ncdf)

# here we set a working directory
setwd("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/working_folder")

###########################################################################################################################################################################
# this is some NetCDF stuff that will aid in the determination of which level to use:
# read the NC
#nc<-open.ncdf("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.cccma_cgcm3_1.20c3m.run1.monthly.hur_a1_20c3m_1_cgcm3.1_t47_1850_2000.nc")
#get the levels
#levels <- as.vector(nc$dim$lev$vals)
###########################################################################################################################################################################

# read in the entire series of data from a netCDF using stack()
ncstack <- stack("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.cccma_cgcm3_1.sresa1b.run1.monthly.hur_a1_sresa1b_1_cgcm3.1_t47_2001_2100.nc") #/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.mpi_echam5.sresa1b.run1.monthly.hur_A1_2001-2050.nc
# read in the 20c3m stack
climstack <- stack("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.cccma_cgcm3_1.20c3m.run1.monthly.hur_a1_20c3m_1_cgcm3.1_t47_1850_2000.nc")

# here we need to be sure that the pclatlong referenced maps have the TRUE lon(0,360)
# -------------------------
# OLD ncstack extent
# class       : Extent
# xmin        : -1.875
# xmax        : 358.125
# ymin        : -89.01354
# ymax        : 89.01354
# -------------------------

# this is used to set the extent to fit a "NORMAL" 0-360 and not the negative to something weird that it was set at.  
e <- extent(c(0,360,-90,90)) 

# here we are taking the newly created extent and applying it to the two GCM netCDF stacks that we just read in
extent(ncstack) <- e
extent(climstack) <- e

# here we are going to use the subset command to grab the files for 1961-1990
#  ** these values only make sense if the 20c3m or historical timeseries dates are monthly 1850-2000
climstack.select <- subset(climstack, 1332:1680)


# loop through the months in the years and calculate the mean monthly climatology for the period of 1961-1990
#  this climatology is created to create anomalies of the GCM futures data using the 20c3m as the baseline period
for(i in 1:12){
	monthList <- seq(i,nlayers(climstack.select),12) # a list of month indexes through subselected series
	if(nchar(i)<2){ month=paste("0",i,sep="")}else{month=paste(i,sep="")} # month naming convention 
	substack <- subset(climstack.select, monthList) # subset the series
	monMean <- mean(substack) # run a mean to get the climatology for the 30 year period
	writeRaster(monMean, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/hur_climatology/","hur_mean_cccma_cgcm3_20c3m_", month,"_1961_1990.tif",sep=""), overwrite=T)
}

print("Completed Mean Monthlies Historical Period")

##########################################################################
# lets do a little RAM cleanup with stuff we no longer require			 #
# rm(climstack) 															 #
# rm(climstack.select)													 # 
# rm(monthList) 															 #
# rm(month) 																 #
# rm(substack) 															 #
# rm(monMean) 															 #
##########################################################################

# list the newly aggregated mean monthly files
l <- list.files("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/hur_climatology/", pattern=".tif$", full.names=T)

# stack up the newly averaged files
climstack.mean <- stack(l)

# list the tiff files from the directory containing the CRU climatologies for the variable being downscaled
#  This line assumes some things about the input data, as in I have already converted these GLOBAL data into
#  Pacific centered Lat long. 
l <- list.files("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/cru_files", pattern="_PCLL.tif", full.names=T)

# use that list object to stack all the files to a stack
crustack <- stack(l)

# read in the extent shapefile for the interpolation analysis
extent.shape <- readShapeSpatial("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/idrisi/iem_downscale_mask_FINAL.shp")

# grab the desired cells from each stack object being used
climstack.mean.desiredCells <- cellsFromExtent(subset(climstack.mean,1,drop=T), extent.shape)
crustack.desiredCells <- cellsFromExtent(subset(crustack,1,drop=T), extent.shape)
ncstack.desiredCells <- cellsFromExtent(subset(ncstack,1,drop=T), extent.shape)
# now extract those cells to new stacks
ncstack.desiredCells.e <- ncstack[ncstack.desiredCells,drop=F]
crustack.desiredCells.e <- crustack[crustack.desiredCells,drop=F]
climstack.mean.desiredCells.e <- climstack.mean[climstack.mean.desiredCells,drop=F]

# this is a list of the years that is used in creating the output naming convention
yearList <- 2001:2100 


# create anomalies 
for(i in 1:12){
	print(paste("	MONTH WORKING = ",i))
	monthList <- seq(i, nlayers(ncstack), 12)
	if(nchar(i)<2){ month=paste("0",i,sep="")}else{month=paste(i,sep="")}
	climstack.current <- subset(climstack.mean.desiredCells.e, i, drop=T)
	
	count=0 # a counter to iterate through the years
	
	for(j in monthList){
		print(paste("anomalies iter ",j))
		count=count+1
		ncstack.current <- subset(ncstack.desiredCells.e, j, drop=T)
		anom <- ncstack.current/climstack.current

		writeRaster(anom, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/anomalies/","hur_mean_cccma_cgcm3_anom_", month,"_",yearList[count],".tif",sep=""), overwrite=T)
	}
}

print("Completed Anomalies Calculation")

# read anomalies back in here and use below in the interpolation
# would be smart here to clean up some of the crazy variable usage if possible.  I feel like there is duplication.

#create a character vector that will store the list of files
l <- character()

# read the anomalies into a stack
# I am doing this is this loop because due to the filename structure they will not "list" chronologically with list.files()
for(i in yearList){
	for(j in 1:12){
		if(nchar(j)<2){ month=paste("0",j,sep="")}else{month=paste(j,sep="")}	

		l <- append(l, paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/anomalies/","hur_mean_cccma_cgcm3_anom_", month,"_",i,".tif",sep=""), after = length(l))
	}
}

ncstack.anom <- stack(l)

in_xy <- coordinates(subset(ncstack.anom,1,drop=T))
out_xy <- coordinates(subset(crustack.desiredCells.e,1,drop=T))

print("Downscaling...")

for(i in 1:12){ # this is a loop that creates an index of the monthly files I want to grab from the ncstack on monthly basis
	monthList <- seq(i, nlayers(ncstack), 12)
	
	if(nchar(i)<2){ month=paste("0",i,sep="")}else{month=paste(i,sep="")}

	print(month)

	cru.current <- subset(crustack.desiredCells.e, i, drop=T)

	count=0 # a little counter

	for(j in monthList){
		count=count+1
		
		ncstack.current <- subset(ncstack.anom, j, drop=T)
		#year.current <- YearList[counter]
		z_in <- getValues(ncstack.current)


		ncstack.anom.spline	<- interp(x=in_xy[,1],y=in_xy[,2],z=z_in,xo=seq(min(out_xy[,1]),max(out_xy[,1]),l=ncol(cru.current)), yo=seq(min(out_xy[,2]),max(out_xy[,2]),l=nrow(cru.current),linear=F))

		# transpose the data here
		nc.interp <- t(ncstack.anom.spline$z)[,nrow(ncstack.anom.spline$z):1]

		# rasterize it

		nc.interp.r <- raster(nc.interp, xmn=xmin(cru.current), xmx=xmax(cru.current), ymn=ymin(cru.current), ymx=ymax(cru.current))


		# write the new raster object to a file
		#tmp.out <- raster(ncstack.anom.spline$z, xmn=xmin(cru.current), xmx=xmax(cru.current), ymn=ymin(cru.current), ymx=ymax(cru.current))

		# multiply the output interpolated raster by the cru 10 min baseline
		downscaled.month <- cru.current*nc.interp.r

		downscaled.month.flip <- cbind(coordinates(downscaled.month), as.numeric(getValues(downscaled.month)))

		if(any(downscaled.month.flip[,1]>180)) downscaled.month.flip[,1][downscaled.month.flip[,1]>180] <- downscaled.month.flip[,1][downscaled.month.flip[,1]>180]-360

		downscaled.month.flip.r <- rasterFromXYZ(downscaled.month.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", digits=8) #downscaled.month.flip.r <- rasterFromXYZ(downscaled.month.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

		downscaled.month.flip.r.v <- getValues(downscaled.month.flip.r)

		ind <- which(downscaled.month.flip.r.v < 0); values(downscaled.month.flip.r)[ind] <- -9999

		writeRaster(downscaled.month.flip.r, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/downscaled/latlong/","hur_mean_cccma_cgcm3_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)

	
		newProj <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +nadgrids=@alaska +datum=NAD83 +units=m +no_defs"

		
		gdalCall <- paste("gdalwarp -t_srs '", newProj, "' -tr 1000 1000 -r bilinear -srcnodata -9999 -dstnodata NoData -co 'COMPRESS=LZW' -cutline '/workspace/UA/malindgren/projects/iem/PHASE2_DATA/Extents/AIEM/AIEM_domain.shp' -overwrite ", paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/downscaled/latlong/","hur_mean_cccma_cgcm3_10min_",month,"_",yearList[count],".tif",sep="")," " , paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/downscaled/AKalbers/","hur_mean_cccma_cgcm3_10min_",month,"_",yearList[count],".tif",sep=""), sep="")

		system(gdalCall)

		#downscaled.month.projected <- raster(paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/downscaled/AKalbers/","hur_mean_cccma_cgcm3_10min_",month,"_",yearList[count],".tif",sep=""))

		#gdalwarp -t_srs "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +nadgrids=@alaska +datum=NAD83 +units=m +no_defs" -tr 1000 1000 -r lanczos -srcnodata -9999 -dstnodata NoData -co "COMPRESS=LZW" /workspace/UA/malindgren/projects/iem/PHASE2_DATA/testing_dir/monthFlip_test.tif /workspace/UA/malindgren/projects/iem/PHASE2_DATA/testing_dir/monthFlip_test_gdal6.tif

		#gdalwarp -t_srs "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +nadgrids=@alaska +datum=NAD83 +units=m +no_defs" -tr 1000 1000 -r lanczos -srcnodata -9999 -dstnodata NoData -co "COMPRESS=LZW" -cutline "/workspace/UA/malindgren/projects/iem/PHASE2_DATA/Extents/AIEM/AIEM_domain.shp" /workspace/UA/malindgren/projects/iem/PHASE2_DATA/testing_dir/monthFlip_test.tif /workspace/UA/malindgren/projects/iem/PHASE2_DATA/testing_dir/monthFlip_test_gdal5.tif




		# 10 min output
		#writeRaster(downscaled.month.flip.r, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/downscaled/AKalbers/","hur_mean_cccma_cgcm3_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)


		#downscaled.month.flip.r.akalb <- projectRaster()

		#######################################################################################################################
		# possibly write in some further interpolation down to 1km using another spline with these outputs...
		# talk to Tom about this...

		#1km.interp <- cbind(coordinates(downscaled.month), getValues(downscaled.month))
		#disaggregate ?

		#resample ?

		#######################################################################################################################
	}
}
	
d.end <- date()

print(d.begin)
print(d.end)



#### TESTING AREA #####

test <- raster(ncstack, layer=1)
test.xyz <- cbind(coordinates(test), getValues(test))
test.xyz
class(test.xyz)
as.data.frame(test.xyz)
test.xyz <- as.data.frame(test.xyz)

colnames(test.xyz) <- c("x","y","z")

m <- matrix(test.xyz$z, nrow=nrow(test), ncol=ncol(test))

if(any(coords[,1]>180)) coords[,1][coords[,1]>180] <- coords[,1][coords[,1]>180]-360

#-----------------------------------------------------------------------------
# here is how we will create the new flipped raster file from an XYZ table

# create an xyz table from 0-360 raster
xyz <- data.frame(coordinates(test),getValues(test))

# set the colnames to something simple
colnames(xyz) <- c("x","y","z")

# turn that 0-360 XYZ table into a -180-180 table
if(any(xyz[,1]>180)) xyz[,1][xyz[,1]>180] <- xyz[,1][xyz[,1]>180]-360

# turn that table into a raster
xyz.r <- rasterFromXYZ(xyz,res=res(test),crs=NA)

#-----------------------------------------------------------------------------

ncstack <- stack("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.cccma_cgcm3_1.sresa1b.run1.monthly.hur_a1_sresa1b_1_cgcm3.1_t47_2001_2100.nc", xmn=0,xmx=360,ymn=-90,ymx=90)

if(any(c.f[,1]>180)) c.f[,1][c.f[,1]>180] <- c.f[,1][c.f[,1]>180]-360


testing <- projectRaster(xxx, newProj)


writeRaster(testing, filename="/workspace/UA/malindgren/projects/iem/PHASE2_DATA/testing_dir/test.tif", overwrite=TRUE)


system("gdalwarp -t_srs '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +nadgrids=@alaska +datum=NAD83 +units=m +no_defs' /workspace/UA/malindgren/projects/iem/PHASE2_DATA/testing_dir/monthFlip_test.tif /workspace/UA/malindgren/projects/iem/PHASE2_DATA/testing_dir/monthFlip_test_gdal.tif")

gdalCall <- paste("gdalwarp -t_srs '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +nadgrids=@alaska +datum=NAD83 +units=m +no_defs' ", <in_file>," " ,<out_file>, sep="")


