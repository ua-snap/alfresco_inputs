d.begin <- date()

# read in the spatial packages needed in this downscaling
require(raster)
require(akima)
require(maptools)
require(rgdal)
require(ncdf)

# here we set a working directory
setwd("/big_storage/malindgren/AIEM/RHUM/working")

outPath <- "/big_storage/malindgren/AIEM/RHUM/"

# read in the extent shapefile for the interpolation analysis
extent.shape <- readShapeSpatial("/big_storage/malindgren/AIEM/RHUM/extent/iem_downscale_mask_FINAL.shp")
model= "cccma_cgcm3_1" # "mpi_echam5"
model.short = unlist(strsplit(model,"_"))[1]

scenario = "a1b"

###########################################################################################################################################################################
# this is some NetCDF stuff that will aid in the determination of which level to use:
# read the NC
#nc<-open.ncdf("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.cccma_cgcm3_1.20c3m.run1.monthly.hur_a1_20c3m_1_cgcm3.1_t47_1850_2000.nc")
#get the levels
#levels <- as.vector(nc$dim$lev$vals)
###########################################################################################################################################################################

# here are a few options that can be set within the raster package
setOptions(tmpdir="/tmp/", progress="text", timer=TRUE, chunksize=1e+24, maxmemory=1e+24)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# this function allows for a very quick stacking of data.  You pass it a list of files from say list.files()
quickStack <- function(f) {
	r <- raster(f[1])
	ln <- extension(basename(f), '')
	s <- stack(r)
	s@layers <- sapply(1:length(f), function(x){ r@file@name = f[x]; r@layernames=ln[x]; r@data@haveminmax=FALSE ; r })
	s@layernames <- ln
	s
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now read in an output mask that will be used in the final step of the processing
out.mask <- raster("/big_storage/malindgren/AIEM/RHUM/mask/ALFRESCO_MASK_1km_AKCanada.tif")

# read in the entire series of data from a netCDF using stack()
hur.future <- stack("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.cccma_cgcm3_1.sresa1b.run1.monthly.hur_a1_sresa1b_1_cgcm3.1_t47_2001_2100.nc") #/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.mpi_echam5.sresa1b.run1.monthly.hur_A1_2001-2050.nc
tas.future <- stack("/Data/Base_Data/Climate/World/GCM_raw/tas/Alaska5/pcmdi.ipcc4.cccma_cgcm3_1.sresa1b.run1.monthly.tas_a1_sresa1b_1_cgcm3.1_t47_2001_2100.nc") #/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.mpi_echam5.sresa1b.run1.monthly.hur_A1_2001-2050.nc

# read in the 20c3m stack
hur.past <- stack("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.cccma_cgcm3_1.20c3m.run1.monthly.hur_a1_20c3m_1_cgcm3.1_t47_1850_2000.nc")
tas.past <- stack("/Data/Base_Data/Climate/World/GCM_raw/tas/Alaska5/pcmdi.ipcc4.cccma_cgcm3_1.20c3m.run1.monthly.tas_a1_20c3m_1_cgcm3.1_t47_1850_2000.nc")

# here we need to be sure that the pclatlong referenced maps have the TRUE lon(0,360)
# -------------------------
# OLD hur.future extent ** when initally imported **
# class       : Extent
# xmin        : -1.875
# xmax        : 358.125
# ymin        : -89.01354
# ymax        : 89.01354
# -------------------------


# this is used to set the extent to fit a "NORMAL" 0-360 and not the negative to something weird that it was set at.  
e <- extent(c(0,360,-90,90)) 

# here we are taking the newly created extent and applying it to the two GCM netCDF stacks that we just read in
extent(hur.future) <- e
extent(hur.past) <- e
extent(tas.future) <- e
extent(tas.past) <- e


#((1961-1850)*12)+1 = 1333  <<<<  this is the begin point on the subset for 1961-1990
#((1990-1850)*12) = 1680

# here we are going to use the subset command to grab the files for 1961-1990
#  ** these values only make sense if the 20c3m or historical timeseries dates are monthly 1850-2000
hur.past.select <- subset(hur.past, 1333:1680)
tas.past.select <- subset(tas.past, 1333:1680)

# loop through the months in the years and calculate the mean monthly climatology for the period of 1961-1990
#  this climatology is created to create anomalies of the GCM futures data using the 20c3m as the baseline period
for(i in 1:12){
	monthList <- seq(i,nlayers(hur.past.select),12) # a list of month indexes through subselected series
	if(nchar(i)<2){ month=paste("0",i,sep="")}else{month=paste(i,sep="")} # month naming convention 
	substack.hur <- subset(hur.past.select, monthList) # subset the series
	substack.tas <- subset(tas.past.select, monthList) # subset the series
	monMean.hur <- mean(substack.hur) # run a mean to get the climatology for the 30 year period
	monMean.tas <- mean(substack.tas) # run a mean to get the climatology for the 30 year period
	writeRaster(monMean.hur, filename=paste(outPath,"climatologies/hur/",model.short,"/","hur_mean_",model,"_20c3m_", month,"_1961_1990.tif",sep=""), overwrite=T)
	writeRaster(monMean.tas, filename=paste(outPath,"climatologies/tas/",model.short,"/","tas_mean_",model,"_20c3m_", month,"_1961_1990.tif",sep=""), overwrite=T)
}

print("Completed Mean Monthlies Historical Period")

# #########################################################################
# lets do a little RAM cleanup with stuff we no longer require			 #
# rm(hur.past) 															 #
# rm(hur.past.select)													 # 
# rm(monthList) 															 #
# rm(month) 																 #
# rm(substack) 															 #
# rm(monMean) 															 #
# #########################################################################

# list the newly aggregated mean monthly files
l <- list.files(paste(outPath,"climatologies/hur/",model.short,"/",sep=""), pattern=".tif$", full.names=T)

# stack up the newly averaged files
hur.past.mean <- stack(l)

l <- list.files(paste(outPath,"climatologies/tas/",model.short,"/",sep=""), pattern=".tif$", full.names=T)

# stack up the newly averaged files
tas.past.mean <- stack(l)

# list the tiff files from the directory containing the CRU climatologies for the variable being downscaled
#  This line assumes some things about the input data, as in I have already converted these GLOBAL data into
#  Pacific centered Lat long. 
# [NEW ML]
# here we read in the CRU TS20 1961-1990 Climatology data as a stack
# list the tiff files from the directory containing the CRU climatologies for the variable being downscaled
#  This line assumes some things about the input data, as in I have already converted these GLOBAL data into
#  Pacific centered Lat long. 
l <- list.files("/big_storage/malindgren/AIEM/RHUM/hur_cru_10min_1961_1990/102212_RASTERIZE/pcll", pattern="_pcll_pcll.tif", full.names=T)

# use that list object to stack all the files to a stack
crustack.hur <- stack(l)

# and the same thing for the temperature data
l <- list.files("/big_storage/malindgren/AIEM/RHUM/tas_cru_10min_1961_1990/cru_10min_splined_extent_pcll", pattern=".tif", full.names=T)
crustack.tas <- stack(l)


# grab the desired cells from each stack object being used
hur.past.mean.desiredCells <- cellsFromExtent(subset(hur.past.mean,1,drop=T), extent.shape)
tas.past.mean.desiredCells <- cellsFromExtent(subset(tas.past.mean,1,drop=T), extent.shape)
crustack.hur.desiredCells <- cellsFromExtent(subset(crustack.hur,1,drop=T), extent.shape)
crustack.tas.desiredCells <- cellsFromExtent(subset(crustack.tas,1,drop=T), extent.shape)
hur.future.desiredCells <- cellsFromExtent(subset(hur.future,1,drop=T), extent.shape)
tas.future.desiredCells <- cellsFromExtent(subset(tas.future,1,drop=T), extent.shape)

# now extract those cells to new stacks
hur.future.desiredCells.e <- hur.future[hur.future.desiredCells,drop=F]
tas.future.desiredCells.e <- tas.future[tas.future.desiredCells,drop=F]
crustack.hur.desiredCells.e <- crustack.hur[crustack.hur.desiredCells,drop=F]
crustack.tas.desiredCells.e <- crustack.tas[crustack.tas.desiredCells,drop=F]
hur.past.mean.desiredCells.e <- hur.past.mean[hur.past.mean.desiredCells,drop=F]
tas.past.mean.desiredCells.e <- tas.past.mean[tas.past.mean.desiredCells,drop=F]


# this is a list of the years that is used in creating the output naming convention
yearList <- 2001:2100 

# create anomalies 
for(i in 1:12){
	print(paste("	MONTH WORKING = ",i))
	monthList <- seq(i, nlayers(hur.future), 12)
	if(nchar(i)<2){ month=paste("0",i,sep="")}else{month=paste(i,sep="")}
	hur.past.current <- subset(hur.past.mean.desiredCells.e, i, drop=T)
	tas.past.current <- subset(tas.past.mean.desiredCells.e, i, drop=T)
	count=0 # a counter to iterate through the years
	
	for(j in monthList){
		print(paste("anomalies iter ",j))
		count=count+1
		hur.future.current <- subset(hur.future.desiredCells.e, j, drop=T)
		tas.future.current <- subset(tas.future.desiredCells.e, j, drop=T)
		anom.hur <- hur.future.current/hur.past.current
		anom.tas <- tas.future.current/tas.past.current

		writeRaster(anom.hur, filename=paste(outPath,"anomalies/hur/",model.short,"/","hur_mean_",model,"_anom_", month,"_",yearList[count],".tif",sep=""), overwrite=T)
		writeRaster(anom.tas, filename=paste(outPath,"anomalies/tas/",model.short,"/","tas_mean_",model,"_anom_", month,"_",yearList[count],".tif",sep=""), overwrite=T)
	}
}

print("Completed Anomalies Calculation")

# read anomalies back in here and use below in the interpolation
# would be smart here to clean up some of the crazy variable usage if possible.  I feel like there is duplication.

#create a character vector that will store the list of files
l.hur <- character()
l.tas <- character()

# read the anomalies into a stack
# I am doing this is this loop because due to the filename structure they will not "list" chronologically with list.files()
for(i in yearList){
	for(j in 1:12){
		if(nchar(j)<2){ month=paste("0",j,sep="")}else{month=paste(j,sep="")}	

		l.hur <- append(l.hur, paste(outPath,"anomalies/hur/",model.short,"/","hur_mean_",model,"_anom_", month,"_",i,".tif",sep=""), after = length(l.hur))
		l.tas <- append(l.tas, paste(outPath,"anomalies/tas/",model.short,"/","tas_mean_",model,"_anom_", month,"_",i,".tif",sep=""), after = length(l.tas))
	}
}

hur.future.anom <- quickStack(l.hur)
tas.future.anom <- quickStack(l.tas)

# ######################################################################################################################################################################################################
# #create a character vector that will store the list of files for the tas data used in the conversion from the relative humidity result to vapor pressure.
# tas.list <- character()

# # read the anomalies into a stack
# # I am doing this is this loop because due to the filename structure they will not "list" chronologically with list.files()
# for(i in yearList){
# 	for(j in 1:12){
# 		if(nchar(j)<2){ month=paste("0",j,sep="")}else{month=paste(j,sep="")}	

# 		tas.list <- append(tas.list, paste("/Data/Base_Data/Climate/AK_CAN_1km_from2km/",model,"/sresa1b/tas/","tas_mean_C_alf_ar4_",model,"_sresa1b_",month,"_",i,".tif",sep=""), after = length(tas.list))
# 	}
# }

# tas.stack <- stack(tas.list)

# ######################################################################################################################################################################################################


in_xy <- coordinates(subset(hur.future.anom,1,drop=T))
out_xy <- coordinates(subset(crustack.hur.desiredCells.e,1,drop=T))

print("Downscaling...")

for(i in 1:12){ # this is a loop that creates an index of the monthly files I want to grab from the hur.future on monthly basis
	monthList <- seq(i, nlayers(hur.future.anom), 12)
	
	if(nchar(i)<2){ month=paste("0",i,sep="")}else{month=paste(i,sep="")}

	print(month)

	cru.hur.current <- subset(crustack.hur.desiredCells.e, i, drop=T)
	cru.tas.current <- subset(crustack.tas.desiredCells.e, i, drop=T)

	count=0 # a little counter

	for(j in monthList){
		count=count+1
		
		hur.future.current <- subset(hur.future.anom, j, drop=T)

		tas.future.current <- subset(tas.future.anom, j, drop=T)

		z_in.hur <- getValues(hur.future.current)
		z_in.tas <- getValues(tas.future.current)

		hur.future.anom.spline	<- interp(x=in_xy[,1],y=in_xy[,2],z=z_in.hur, xo=seq(min(out_xy[,1]),max(out_xy[,1]), l=ncol(cru.hur.current)), yo=seq(min(out_xy[,2]),max(out_xy[,2]), l=nrow(cru.hur.current)))
		
		tas.future.anom.spline	<- interp(x=in_xy[,1],y=in_xy[,2],z=z_in.tas, xo=seq(min(out_xy[,1]),max(out_xy[,1]), l=ncol(cru.tas.current)), yo=seq(min(out_xy[,2]),max(out_xy[,2]), l=nrow(cru.tas.current)))
		
		# transpose the data here
		mat.hur <- matrix(hur.future.anom.spline$z,nrow=nrow(cru.hur.current),ncol=ncol(cru.hur.current),byrow=T)[nrow(cru.hur.current):1,]
		mat.tas <- matrix(tas.future.anom.spline$z,nrow=nrow(cru.tas.current),ncol=ncol(cru.tas.current),byrow=T)[nrow(cru.tas.current):1,]

		# nc.interp.hur <- t(hur.future.anom.spline$z)[,nrow(hur.future.anom.spline$z):1]
		# nc.interp.tas <- t(tas.future.anom.spline$z)[,nrow(tas.future.anom.spline$z):1]

		# rasterize it
		nc.interp.hur.r <- raster(mat.hur, xmn=xmin(cru.hur.current), xmx=xmax(cru.hur.current), ymn=ymin(cru.hur.current), ymx=ymax(cru.hur.current))
		nc.interp.tas.r <- raster(mat.tas, xmn=xmin(cru.tas.current), xmx=xmax(cru.tas.current), ymn=ymin(cru.tas.current), ymx=ymax(cru.tas.current))

		# write the new raster object to a file
		#tmp.out <- raster(hur.future.anom.spline$z, xmn=xmin(cru.hur.current), xmx=xmax(cru.hur.current), ymn=ymin(cru.hur.current), ymx=ymax(cru.hur.current))

		# multiply the output interpolated raster by the cru 10 min baseline
		downscaled.month.hur <- cru.hur.current*nc.interp.hur.r
		downscaled.month.tas <- cru.tas.current*nc.interp.tas.r

		# here we create an xyz table which will be used to flip the coordinates back to greenwich centered Latlong
		downscaled.month.hur.flip <- cbind(coordinates(downscaled.month.hur), as.numeric(getValues(downscaled.month.hur)))
		downscaled.month.tas.flip <- cbind(coordinates(downscaled.month.tas), as.numeric(getValues(downscaled.month.tas)))

		# this line changes the values of longitude from 0-360 to -180 - 180
		if(any(downscaled.month.hur.flip[,1]>180)) downscaled.month.hur.flip[,1][downscaled.month.hur.flip[,1]>180] <- downscaled.month.hur.flip[,1][downscaled.month.hur.flip[,1]>180]-360
		if(any(downscaled.month.tas.flip[,1]>180)) downscaled.month.tas.flip[,1][downscaled.month.tas.flip[,1]>180] <- downscaled.month.tas.flip[,1][downscaled.month.tas.flip[,1]>180]-360

		# turn that xyz table into a raster again with the corrected projection system
		downscaled.month.hur.flip.r <- rasterFromXYZ(downscaled.month.hur.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", digits=8) #downscaled.month.hur.flip.r <- rasterFromXYZ(downscaled.month.hur.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
		downscaled.month.tas.flip.r <- rasterFromXYZ(downscaled.month.tas.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", digits=8) #downscaled.month.tas.flip.r <- rasterFromXYZ(downscaled.month.tas.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

		# here we get the values of the current downscaled layer
		downscaled.month.hur.flip.r.v <- getValues(downscaled.month.hur.flip.r)
		
		# now we ask some questions about the values in the map.  The acceptable range is 0-100 ** although there are instances where there is >100 
		# here we are asking it to turn all of these out of bounds values to a NoData value of -9999
		ind <- which(downscaled.month.hur.flip.r.v < 0); values(downscaled.month.hur.flip.r)[ind] <- -9999

		# convert back to vapor pressure
		esa = 6.112*exp(17.62*downscaled.month.tas.flip.r/(243.12+downscaled.month.tas.flip.r))
		vapor.ts31.10min = (downscaled.month.hur.flip.r*esa)/100

		values(downscaled.month.hur.flip.r)[which(getValues(downscaled.month.hur.flip.r) > 100)] <- 95

		# here we write out the greenwich-centered version of the map
		writeRaster(trim(downscaled.month.hur.flip.r), filename=paste(outPath,"downscaled/",scenario,"/hur/",model.short,"/","hur_mean_",model,"_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)
		writeRaster(trim(downscaled.month.tas.flip.r), filename=paste(outPath,"downscaled/",scenario,"/tas/",model.short,"/","tas_mean_",model,"_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)
		writeRaster(trim(vapor.ts31.10min), filename=paste(outPath,"downscaled/",scenario,"/vap/",model.short,"/","vap_mean_",model,"_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)

		# this is the proj4 string to reproject into Alaska Albers 
		newProj <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +nadgrids=@alaska +datum=NAD83 +units=m +no_defs"

		# In this instance I have found that the Raster package has a dummed down version of the proj4 lib and it is giving strange results
		# therefore I am using a system call from R to the GDAL library to do the reproject/resample/clipping of the data...
		# gdalCall <- paste("gdalwarp -t_srs '", newProj, "' -tr 1000 1000 -r bilinear -srcnodata -9999 -dstnodata NoData -co 'COMPRESS=LZW' -cutline '/big_storage/malindgren/AIEM/RHUM/extent/AIEM_domain.shp' -overwrite ", paste(outPath,"downscaled/",scenario,"/",model.short,"/latlong/","hur_mean_",model,"_10min_",month,"_",yearList[count],".tif",sep="")," " ,paste(outPath,"downscaled/",scenario,"/",model.short,"/akalbers/","hur_mean_",model,"_10min_",month,"_",yearList[count],".tif",sep=""), sep="")
		gdalCall <- paste("gdalwarp -t_srs '", newProj, "' -tr 1000 1000 -r bilinear -srcnodata -9999 -dstnodata NoData -co 'COMPRESS=LZW' -overwrite ", paste(outPath,"downscaled/",scenario,"/hur/",model.short,"/10min/","hur_mean_",model,"_10min_",month,"_",yearList[count],".tif",sep="")," " ,paste(outPath,"downscaled/",scenario,"/hur/",model.short,"/1km/","hur_mean_",model,"_1km_",month,"_",yearList[count],".tif",sep=""), sep="")

		# this is how we issue a system call in R.
		system(gdalCall)

		gdalCall <- paste("gdalwarp -t_srs '", newProj, "' -tr 1000 1000 -r bilinear -srcnodata -9999 -dstnodata NoData -co 'COMPRESS=LZW' -overwrite ", paste(outPath,"downscaled/",scenario,"/tas/",model.short,"/10min/","tas_mean_",model,"_10min_",month,"_",yearList[count],".tif",sep="")," " ,paste(outPath,"downscaled/",scenario,"/tas/",model.short,"/1km/","tas_mean_",model,"_1km_",month,"_",yearList[count],".tif",sep=""), sep="")

		system(gdalCall)

		gdalCall <- paste("gdalwarp -t_srs '", newProj, "' -tr 1000 1000 -r bilinear -srcnodata -9999 -dstnodata NoData -co 'COMPRESS=LZW' -overwrite ", paste(outPath,"downscaled/",scenario,"/vap/",model.short,"/10min/","vap_mean_",model,"_10min_",month,"_",yearList[count],".tif",sep="")," " ,paste(outPath,"downscaled/",scenario,"/vap/",model.short,"/1km/","vap_mean_",model,"_1km_",month,"_",yearList[count],".tif",sep=""), sep="")

		system(gdalCall)

		# now we need to mask them
		# read back in the 1km hur
		hur.mask <- raster(paste(outPath,"downscaled/",scenario,"/hur/",model.short,"/1km/","hur_mean_",model,"_1km_",month,"_",yearList[count],".tif",sep=""))
		# read back in the 1km vapor
		vap.mask <- raster(paste(outPath,"downscaled/",scenario,"/tas/",model.short,"/1km/","tas_mean_",model,"_1km_",month,"_",yearList[count],".tif",sep=""))
		# read back in the 1km tas
		tas.mask <- raster(paste(outPath,"downscaled/",scenario,"/vap/",model.short,"/1km/","vap_mean_",model,"_1km_",month,"_",yearList[count],".tif",sep=""))

		# and now we mask them and write 'em out
		mask(hur.mask, out.mask, filename=paste(outPath,"downscaled/",scenario,"/hur/",model.short,"/1km_mask/","hur_mean_",model,"_1km_",month,"_",yearList[count],".tif",sep=""), options="COMPRESS=LZW")
		mask(vap.mask, out.mask, filename=paste(outPath,"downscaled/",scenario,"/tas/",model.short,"/1km_mask/","tas_mean_",model,"_1km_",month,"_",yearList[count],".tif",sep=""), options="COMPRESS=LZW")
		mask(tas.mask, out.mask, filename=paste(outPath,"downscaled/",scenario,"/vap/",model.short,"/1km_mask/","vap_mean_",model,"_1km_",month,"_",yearList[count],".tif",sep=""), options="COMPRESS=LZW")
	}
}

d.end <- date()

print(d.begin)
print(d.end)



