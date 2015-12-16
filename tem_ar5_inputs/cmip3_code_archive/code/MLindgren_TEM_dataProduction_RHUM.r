d.begin <- date()

# read in the spatial packages needed in this downscaling
require(raster)
require(akima)
require(maptools)
require(rgdal)
require(ncdf)

# here we set a working directory
setwd("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/working_folder")

outPath <- "/big_storage/malindgren/AIEM/RHUM/"

# read in the extent shapefile for the interpolation analysis
extent.shape <- readShapeSpatial("/big_storage/malindgren/AIEM/RHUM/extent/iem_downscale_mask_FINAL.shp")
model= "cccma_cgcm3_1" # "mpi_echam5"
model.short = unlist(strsplit(model,"_"))[1]

scenario = "20c3m"
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
# read in the needed data 


# here we need to be sure that the pclatlong referenced maps have the TRUE lon(0,360)
# -------------------------
# OLD ncstack extent ** when initally imported **
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
	writeRaster(monMean, filename=paste(outPath,"climatologies/",scenario,"/",model.short,"/","hur_mean_",model,"_20c3m_", month,"_1961_1990.tif",sep=""), overwrite=T)
}

print("Completed Mean Monthlies Historical Period")

# #########################################################################
# lets do a little RAM cleanup with stuff we no longer require			 #
# rm(climstack) 															 #
# rm(climstack.select)													 # 
# rm(monthList) 															 #
# rm(month) 																 #
# rm(substack) 															 #
# rm(monMean) 															 #
# #########################################################################

# list the newly aggregated mean monthly files
l <- list.files(paste(outPath,"climatologies/",scenario,"/",model.short,"/",sep=""), pattern=".tif$", full.names=T)

# stack up the newly averaged files
climstack.mean <- stack(l)

# list the tiff files from the directory containing the CRU climatologies for the variable being downscaled
#  This line assumes some things about the input data, as in I have already converted these GLOBAL data into
#  Pacific centered Lat long. 
# [NEW ML]
# here we read in the CRU TS20 1961-1990 Climatology data as a stack
# this data have already been pre-processed into (located on the data drive at snap
crustack <- stack(list.files("/Data/Base_Data/Climate/Canada/CRU_Climatology/Relative_Humidity/raster_grids_pacificCenteredGlobal", pattern=".tif", full.names=T))


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

		writeRaster(anom, filename=paste(outPath,"anomalies","/",scenario,"/",model.short,"/","hur_mean_",model,"_anom_", month,"_",yearList[count],".tif",sep=""), overwrite=T)
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

		l <- append(l, paste(outPath,"anomalies","/",scenario,"/",model.short,"/","hur_mean_",model,"_anom_", month,"_",i,".tif",sep=""), after = length(l))
	}
}

ncstack.anom <- stack(l)

######################################################################################################################################################################################################
#create a character vector that will store the list of files for the tas data used in the conversion from the relative humidity result to vapor pressure.
tas.list <- character()

# read the anomalies into a stack
# I am doing this is this loop because due to the filename structure they will not "list" chronologically with list.files()
for(i in yearList){
	for(j in 1:12){
		if(nchar(j)<2){ month=paste("0",j,sep="")}else{month=paste(j,sep="")}	

		tas.list <- append(tas.list, paste("/Data/Base_Data/Climate/AK_CAN_1km_from2km/",model,"/sresa1b/tas/","tas_mean_C_alf_ar4_",model,"_sresa1b_",month,"_",i,".tif",sep=""), after = length(tas.list))
	}
}

tas.stack <- stack(tas.list)

######################################################################################################################################################################################################


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

		# here we subset the stack of the tas data for 
		tas.stack.current <- subset(tas.stack, j, drop=T)

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

		# here we create an xyz table which will be used to flip the coordinates back to greenwich centered Latlong
		downscaled.month.flip <- cbind(coordinates(downscaled.month), as.numeric(getValues(downscaled.month)))

		# this line changes the values of longitude from 0-360 to -180 - 180
		if(any(downscaled.month.flip[,1]>180)) downscaled.month.flip[,1][downscaled.month.flip[,1]>180] <- downscaled.month.flip[,1][downscaled.month.flip[,1]>180]-360

		# turn that xyz table into a raster again with the corrected projection system
		downscaled.month.flip.r <- rasterFromXYZ(downscaled.month.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", digits=8) #downscaled.month.flip.r <- rasterFromXYZ(downscaled.month.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

		# here we get the values of the current downscaled layer
		downscaled.month.flip.r.v <- getValues(downscaled.month.flip.r)

		# now we ask some questions about the values in the map.  The acceptable range is 0-100 ** although there are instances where there is >100 
		# here we are asking it to turn all of these out of bounds values to a NoData value of -9999
		ind <- which(downscaled.month.flip.r.v < 0); values(downscaled.month.flip.r)[ind] <- -9999

		# here we write out the greenwich-centered version of the map
		writeRaster(downscaled.month.flip.r, filename=paste(outPath,"downscaled/",scenario,"/",model.short,"/latlong/","hur_mean_",model,"_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)

		# this is the proj4 string to reproject into Alaska Albers 
		newProj <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +nadgrids=@alaska +datum=NAD83 +units=m +no_defs"

		# In this instance I have found that the Raster package has a dummed down version of the proj4 lib and it is giving strange results
		# therefore I am using a system call from R to the GDAL library to do the reproject/resample/clipping of the data...
		gdalCall 	<- paste("gdalwarp -t_srs '", newProj, "' -tr 1000 1000 -r bilinear -srcnodata -9999 -dstnodata NoData -co 'COMPRESS=LZW' -cutline '/big_storage/malindgren/AIEM/RHUM/extent/AIEM_domain.shp' -overwrite ", paste(outPath,"downscaled/",scenario,"/",model.short,"/latlong/","hur_mean_",model,"_10min_",month,"_",yearList[count],".tif",sep="")," " ,paste(outPath,"downscaled/",scenario,"/",model.short,"/akalbers/","hur_mean_",model,"_10min_",month,"_",yearList[count],".tif",sep=""), sep="")

		# this is how we issue a system call in R.
		system(gdalCall)

		rhum.downscaled <- raster(paste(outPath,"downscaled/",scenario,"/",model.short,"/akalbers/","hur_mean_",model,"_10min_",month,"_",yearList[count],".tif",sep=""))

		#######################################################################################################################
		# this is the conversion as I got it from Fengming and the TEM group many moons ago
		# esa = 6.112 exp(17.62*Tair/(243.12 + Tair)),  where Tair in °C and >0

		# Conversion Procedure from Fengming vetted by Steph McAfee
		# source: WMO_CIMO_Guide-7th_Edition-2008
		# this equation was modified to isolate the EA by Michael Lindgren (malindgren@alaska.edu)
		ESA = 6.112 ^(22.46*tas.stack.current / (272.62 + tas.stack.current))
		EA = (ESA*rhum.downscaled)/100

		# write out the new raster file
		writeRaster(EA, filename=paste(outPath,"downscaled/vapor/",scenario,"/",model.short,"/","vap_mean_",model,"_1km_",month,"_",yearList[count],".tif",sep=""), overwrite=T)
		
		#######################################################################################################################
	}
}

d.end <- date()

print(d.begin)
print(d.end)

