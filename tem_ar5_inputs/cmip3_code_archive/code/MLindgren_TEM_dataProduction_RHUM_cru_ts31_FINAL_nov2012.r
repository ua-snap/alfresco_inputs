d.begin <- date()

# read in the spatial packages needed in this downscaling
require(raster)
require(akima)
require(maptools)
require(rgdal)
require(ncdf)

# here we set a working directory
setwd("/big_storage/malindgren/AIEM/RHUM/working")

# this is the base output path that is the container for all of the outputs
outPath <- "/big_storage/malindgren/AIEM/RHUM/"

# here are a few options that can be set within the raster package
setOptions(tmpdir="/tmp/", progress="text", timer=TRUE, chunksize=1e+24, maxmemory=1e+24)

# here we bring in some functions to be used in the processing...
source("/big_storage/malindgren/AIEM/functions/CalcAnom.r")


##########################################################################################################################################################
# this is some NetCDF stuff that will aid in the determination of which level to use:
# these lines can be uncommented and run with any NetCDF file that a user may need some metadata about.  
# Namely the levels and varnames within the .nc
# read the NC
# nc<-open.ncdf("<path_to_netcdf.nc>")
# get the levels
# levels <- as.vector(nc$dim$lev$vals)
##########################################################################################################################################################

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

# read in the entire series of data from a netCDF using brick() because the brick function allows for fast reading in of lots of files
vap.stack <- brick("/Data/Base_Data/Climate/World/CRU_grids/CRU_TS31/cru_ts_3_10.1901.2009.vap.dat.nc")
tas.stack <- brick("/Data/Base_Data/Climate/World/CRU_grids/CRU_TS31/cru_ts_3_10.1901.2009.tmp.dat.nc")

# list the tiff files from the directory containing the CRU climatologies for the variable being downscaled
#  This line assumes some things about the input data, as in I have already converted these GLOBAL data into
#  Pacific centered Lat long. 
# l <- list.files("/big_storage/malindgren/AIEM/RHUM/hur_cru_10min_1961_1990/pcll_spline_final", pattern=".tif", full.names=T)

# use that list object to stack all the files to a stack
crustack.hur <- stack(list.files("/big_storage/malindgren/AIEM/RHUM/hur_cru_10min_1961_1990/pcll_spline_final", pattern=".tif", full.names=T))

# l <- list.files("/big_storage/malindgren/AIEM/RHUM/tas_cru_10min_1961_1990/pcll_spline_final", pattern=".tif", full.names=T)
crustack.tas <- stack(list.files("/big_storage/malindgren/AIEM/RHUM/tas_cru_10min_1961_1990/pcll_spline_final", pattern=".tif", full.names=T))

# read in the extent shapefile for the interpolation analysis (in Pacific Centered LatLong)
extent.shape <- extent(c(160, 310, 30, 84)) # OLD >>> extent(c(166.6012,305.4721,38.023,83.4692)) #readShapeSpatial("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/WORKING_TEM/Extents/Extents_TEM_dataCreate/iem_downscale_mask_FINAL.shp")

# this is a list of the years that is used in creating the output naming convention
yearList <- 1901:2009 

# THIS CONVERSION IS HAPPENING HERE BASED ON SUGGESTION FROM STEPH MCAFEE.  IT IS PARTICULARLY IMPORTANT TO HAVE IT HERE 
#  SINCE THERE IS A REPROJECT OF THE STACK OBJECT TO PACIFIC CENTERED LATLONG AND BETTER NOT TO REPROJECT MORE THAN IS NECESSARY.
############################ CONVERSION TO RHUM ##########################################################
# this is the conversion as I got it from Fengming and the TEM group many moons ago
# TEM actually uses vapor pressure to drive the model. So basically the rel hum I 
# sent to you was converted from vapor pressure. The conversion is conducted as following:
# (1) calculate saturated vapor pressure (esa, hPa or mbar) according to Meteorological Instruments and 
# 	  Methods of Observation (CIMO Guide) (WMO, 2008):
#         esa = 6.112 exp(17.62*Tair/(243.12 + Tair)),  where Tair in °C and >0; 
#        OR, 
#         esa = 6.112 exp(22.46 *Tair/(272.62 + Tair)), where Tair  in °C and <=0

# (2) relative humility RH:
#        RH = ea/esa * 100, where ea is the vapor pressure in hPa or mbar.

# Conversion Procedure from Fengming vetted by Steph McAfee
# source: WMO_CIMO_Guide-7th_Edition-2008
# this equation was modified to isolate the EA by Michael Lindgren (malindgren@alaska.edu)
# ESA = 6.112 ^((22.46*tas.stack.desiredCells.e) / (272.62 + tas.stack.desiredCells.e))
# hur.stack = vap.stack.desiredCells.e/ESA*100
print("*****   CONVERTING VAPOR PRESSURE INTO RELATIVE HUMIDITY   ************")

# back to the old conversion
# esa.stack = 6.112*exp(17.62*tas.stack/(243.12+tas.stack))
# hur.stack <- vap.stack/esa.stack * 100

# COMMENTED OUT SINCE ALREADY CREATED

# get the values in a matrix of all of the stacked data
# tas.mat <- getValues(tas.stack)
# vap.mat <- getValues(vap.stack)

# # a list of names from the rasterStack object
# tas.stack.names <- tas.stack@layernames
# tas.stack.names <- strsplit(tas.stack.names, ".", fixed=T)

# # this is used to flip the rasters in this loop as well
# vars.xy <- coordinates(tas.stack)

# # now lets make those coordinates pclatlong 
# if(any(vars.xy[,1] < 0))  vars.xy[vars.xy[,1]<0,1] <- vars.xy[vars.xy[,1]<0,1]+360

# # empty matrices to fill with the calculated output
# esa.mat <- matrix(NA, nrow=nrow(tas.mat), ncol=ncol(tas.mat))
# hur.mat <- matrix(NA, nrow=nrow(tas.mat), ncol=ncol(tas.mat))


# # perform the needed calculation to convert to hur
# esa.mat <- 6.112*exp(17.62*tas.mat/(243.12+tas.mat))
# hur.mat <- vap.mat/esa.mat * 100

# # lets re-order the esa data here using the longitudes 0-360
# hur.mat <- cbind(vars.xy,hur.mat)
# # now we will sort that flipped data
# hur.mat <- hur.mat[order(hur.mat[,1]),]

# hur.mat <- hur.mat[,3:ncol(hur.mat)]

# hur.stack <- brick(nl=nlayers(vap.stack),nrows=nrow(vap.stack),ncols=ncol(vap.stack),xmn=0,xmx=360,ymn=ymin(vap.stack),ymx=ymax(vap.stack),crs=projection(vap.stack))

# hur.stack <- setValues(hur.stack, hur.mat)

### END COMMENTED OUT SINCE ALREADY CREATED

# &&&&&&&&&&&&   IF ABOVE IS COMMENTED THEN RUN THIS PART   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# use these as iterators to create the files you want
monList <- 1:12
yearList <- 1901:2009

hur.l <- character()
tas.l <- character()

# we are using this to read in the data if it is commented out above
iter=0
for(y in yearList){
	for(mon in monList){
		
		if(nchar(mon)<2){ mon=paste("0",mon,sep="")}else{mon=paste(mon,sep="")}	
		
		hur.l <- append(hur.l, paste(outPath,"cru_historical_1901_2009/pacific/hur/","hur_total_mm_cru_TS31_pcLL_", mon,"_",y,".tif",sep=""), after = length(hur.l))
		tas.l <- append(tas.l, paste(outPath,"cru_historical_1901_2009/pacific/tas/","tas_total_mm_cru_TS31_pcLL_", mon,"_",y,".tif",sep=""), after = length(tas.l))
	}
}

# this may no longer be needed <<< CHECK IT AFTER COMPLETION
# this line will re-initialize vap.stack with the newly flipped values
hur.stack <- quickStack(hur.l)
tas.stack <- quickStack(tas.l)
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

print("*****   SELECTING RASTER CELLS USED IN ANALYSIS   ************")
# Here we extract the cells of interest for this calculation from the 2 bricks 
hur.stack.desiredCells <- cellsFromExtent(subset(hur.stack,1,drop=T), extent.shape)
tas.stack.desiredCells <- cellsFromExtent(subset(tas.stack,1,drop=T), extent.shape)
crustack.hur.desiredCells <- cellsFromExtent(subset(crustack.hur,1,drop=T), extent.shape)
crustack.tas.desiredCells <- cellsFromExtent(subset(crustack.tas,1,drop=T), extent.shape)

hur.stack.desiredCells.e <- hur.stack[hur.stack.desiredCells,drop=F]
tas.stack.desiredCells.e <- tas.stack[tas.stack.desiredCells,drop=F]
crustack.hur.desiredCells.e <- crustack.hur[crustack.hur.desiredCells,drop=F]
crustack.tas.desiredCells.e <- crustack.tas[crustack.tas.desiredCells,drop=F]


print("*****   CREATING MONTHLY CLIMATOLOGY & CALCULATING ANOMALIES   ************************")

# here we use the calcAnom() function I wrote and sourced in above
hur.stack.anom <- calcAnom(hur.stack.desiredCells.e,721,1080,absolute=FALSE)
tas.stack.anom <- calcAnom(tas.stack.desiredCells.e,721,1080,absolute=TRUE)


# some grids we will use in the interpolation step
in_xy <- coordinates(subset(hur.stack.anom,1,drop=T))
out_xy <- coordinates(subset(crustack.hur.desiredCells.e,1,drop=T))
# for the 1km interpolation
out_xy_1km <- coordinates(out.mask)

print("*****    PERFORMING DOWNSCALING OF TIMESERIES   ******************")

for(i in 1:12){ # this is a loop that creates an index of the monthly files I want to grab from the vap.stack on monthly basis
	monthList <- seq(i, nlayers(hur.stack.anom), 12)
	
	if(nchar(i)<2){ month=paste("0",i,sep="")}else{month=paste(i,sep="")}

	print(month)

	cru.hur.current <- subset(crustack.hur.desiredCells.e, i, drop=T)
	cru.tas.current <- subset(crustack.tas.desiredCells.e, i, drop=T)

	count=0 # a little counter

	for(j in monthList){
		count=count+1
	
		# downscale the relative humidity	
		hur.stack.current <- subset(hur.stack.anom, j, drop=T)
		in_xyz <- na.omit(data.frame(in_xy, getValues(hur.stack.current)))

		# here we spline the data
		hur.stack.anom.spline <- interp(x=in_xyz[,1],y=in_xyz[,2],z=in_xyz[,3],xo=seq(min(out_xy[,1]),max(out_xy[,1]),l=ncol(cru.hur.current)), yo=seq(min(out_xy[,2]),max(out_xy[,2]),l=nrow(cru.hur.current)))

		# transpose the data here 
		nc.interp <- matrix(hur.stack.anom.spline$z,nrow=nrow(cru.hur.current),ncol=ncol(cru.hur.current),byrow=T)[nrow(cru.hur.current):1,]
		
		#nc.interp <- t(hur.stack.anom.spline$z)[,nrow(hur.stack.anom.spline$z):1]

		# rasterize it
		nc.interp.r <- raster(nc.interp, xmn=xmin(cru.hur.current), xmx=xmax(cru.hur.current), ymn=ymin(cru.hur.current), ymx=ymax(cru.hur.current))

		# write the new raster object to a file
		#tmp.out <- raster(vap.stack.anom.spline$z, xmn=xmin(cru.current), xmx=xmax(cru.current), ymn=ymin(cru.current), ymx=ymax(cru.current))

		# multiply the output interpolated raster by the cru 10 min baseline
		downscaled.month.hur <- cru.hur.current*nc.interp.r

		# # here we create an xyz table which will be used to flip the coordinates back to greenwich centered Latlong
		downscaled.month.flip <- cbind(coordinates(downscaled.month.hur), as.numeric(getValues(downscaled.month.hur)))
		
		# # this line changes the values of longitude from 0-360 to -180 - 180
		if(any(downscaled.month.flip[,1]>180)) downscaled.month.flip[,1][downscaled.month.flip[,1]>180] <- downscaled.month.flip[,1][downscaled.month.flip[,1]>180]-360

		# # turn that xyz table into a raster again with the corrected projection system
		downscaled.month.flip.hur.r <- rasterFromXYZ(downscaled.month.flip,res=res(crustack.hur.desiredCells.e),crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", digits=0) #downscaled.month.flip.r <- rasterFromXYZ(downscaled.month.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")


		### downscale the temperature
		tas.stack.current <- subset(tas.stack.anom, j, drop=T)

		in_xyz <- na.omit(data.frame(in_xy, getValues(tas.stack.current)))

		tas.stack.anom.spline <- interp(x=in_xyz[,1],y=in_xyz[,2],z=in_xyz[,3],xo=seq(min(out_xy[,1]),max(out_xy[,1]),l=ncol(cru.tas.current)), yo=seq(min(out_xy[,2]),max(out_xy[,2]),l=nrow(cru.tas.current)))

		# transpose the data here 
		# check this result!
		nc.interp <- matrix(hur.stack.anom.spline$z,nrow=nrow(cru.hur.current),ncol=ncol(cru.hur.current),byrow=T)[nrow(cru.hur.current):1,]

		#nc.interp <- t(tas.stack.anom.spline$z)[,nrow(tas.stack.anom.spline$z):1]

		# rasterize it
		nc.interp.r <- raster(nc.interp, xmn=xmin(cru.tas.current), xmx=xmax(cru.tas.current), ymn=ymin(cru.tas.current), ymx=ymax(cru.tas.current))

		# write the new raster object to a file
		#tmp.out <- raster(vap.stack.anom.spline$z, xmn=xmin(cru.current), xmx=xmax(cru.current), ymn=ymin(cru.current), ymx=ymax(cru.current))

		# multiply the output interpolated raster by the cru 10 min baseline
		downscaled.month.tas <- cru.tas.current+nc.interp.r


		# # here we create an xyz table which will be used to flip the coordinates back to greenwich centered Latlong
		downscaled.month.flip <- cbind(coordinates(downscaled.month.tas), as.numeric(getValues(downscaled.month.tas)))
		
		# # this line changes the values of longitude from 0-360 to -180 - 180
		if(any(downscaled.month.flip[,1]>180)) downscaled.month.flip[,1][downscaled.month.flip[,1]>180] <- downscaled.month.flip[,1][downscaled.month.flip[,1]>180]-360

		# # turn that xyz table into a raster again with the corrected projection system
		downscaled.month.flip.tas.r <- rasterFromXYZ(downscaled.month.flip, res=res(crustack.hur.desiredCells.e), crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", digits=0) #downscaled.month.flip.r <- rasterFromXYZ(downscaled.month.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

		#### END TAS

		# # here we get the values of the current downscaled layer
		downscaled.month.flip.hur.r.v <- getValues(downscaled.month.flip.hur.r)
		downscaled.month.flip.tas.r.v <- getValues(downscaled.month.flip.tas.r)

		# # now we ask some questions about the values in the map.  The acceptable range is 0-100 ** although there are instances where there is >100 
		# # here we are asking it to turn all of these out of bounds values to a NoData value of -9999
		# ind <- which(downscaled.month.flip.hur.r.v < 0); values(downscaled.month.flip.hur.r)[ind] <- -9999
		# ind <- which(downscaled.month.flip.tas.r.v < 0); values(downscaled.month.flip.tas.r)[ind] <- -9999

		# now we change the >100 values to 95 as per steph's okaying it
		values(downscaled.month.flip.hur.r)[which(values(downscaled.month.flip.hur.r) > 100)] <- 95

		print("*****  COVNERTING RELATIVE HUMIDITY BACK TO VAPOR PRESSURE   *****************")

		# convert back to vapor pressure
		esa = 6.112*exp(17.62*downscaled.month.flip.tas.r/(243.12+downscaled.month.flip.tas.r))
		vapor.ts31.10min = (downscaled.month.flip.hur.r*esa)/100

		# lets trim these data
		downscaled.month.flip.hur.r <- trim(downscaled.month.flip.hur.r)
		downscaled.month.flip.tas.r <- trim(downscaled.month.flip.tas.r)
		vapor.ts31.10min <- trim(vapor.ts31.10min)

		# here we write out the greenwich-centered version of the map at the 10 minute resolution
		writeRaster(downscaled.month.flip.hur.r, filename=paste(outPath,"downscaled/hur_cru_ts31/10min/","hur_cru_ts31_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)
		writeRaster(downscaled.month.flip.tas.r, filename=paste(outPath,"downscaled/tas_cru_ts31/10min/","tas_cru_ts31_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)
		writeRaster(vapor.ts31.10min, filename=paste(outPath,"downscaled/vap_cru_ts31/10min/","vap_cru_ts31_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)


		# this is the proj4 string to reproject into Alaska Albers 
		akalb <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +nadgrids=@alaska +datum=NAD83 +units=m +no_defs"

		# here we reproject the extents of the rasters to alaska albers 
		hur.akalb.ext <- projectExtent(downscaled.month.flip.hur.r, akalb)
		tas.akalb.ext <- projectExtent(downscaled.month.flip.tas.r, akalb)
		vap.akalb.ext <- projectExtent(vapor.ts31.10min, akalb)

		# then we take those extents and reproject the rasters themselves
		hur.akalb <- projectRaster(downscaled.month.flip.hur.r,hur.akalb.ext, method='bilinear')
		tas.akalb <- projectRaster(downscaled.month.flip.tas.r,tas.akalb.ext, method='bilinear')
		vap.akalb <- projectRaster(vapor.ts31.10min,vap.akalb.ext, method='bilinear')

		akalb.vars <- c(hur.akalb,tas.akalb,vap.akalb)
		akalb.vars.names <- c(paste(outPath,"downscaled/hur_cru_ts31/1km_mask/","hur_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""),paste(outPath,"downscaled/vap_cru_ts31/1km_mask/","vap_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""),paste(outPath,"downscaled/tas_cru_ts31/1km_mask/","tas_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""))
		
		for(n in 1:length(akalb.vars)){
			print(paste("working on: ", akalb.vars.names[n], sep=""))
			# it may be the best practice here to just interpolate the data to the desired resolution so that we can remove some of hte data issues inherent in coming from 0.5 to 10 min
			in_xyz <- na.omit(data.frame(coordinates(akalb.vars[n][[1]]), getValues(akalb.vars[n][[1]])))
			akalb.1km.spline <- interp(x=in_xyz[,1],y=in_xyz[,2],z=in_xyz[,3],xo=seq(min(out_xy_1km[,1]),max(out_xy_1km[,1]),l=ncol(out.mask)), yo=seq(min(out_xy_1km[,2]),max(out_xy_1km[,2]),l=nrow(out.mask)))
			# transpose the data here 
			nc.interp <- matrix(akalb.1km.spline$z,nrow=nrow(out.mask),ncol=ncol(out.mask),byrow=T)[nrow(out.mask):1,]
			# rasterize it
			nc.interp.r <- raster(nc.interp, xmn=xmin(out.mask), xmx=xmax(out.mask), ymn=ymin(out.mask), ymx=ymax(out.mask))
			
			print("mask")
			# crop the data to the extent of the output mask
			mask.c <- crop(nc.interp.r,out.mask)
			# then set the damn extent to the out.mask
			extent(mask.c) <- extent(out.mask)
			# now lets mask the output to match ALFRESCO stuff
			mask(mask.c, out.mask, filename=akalb.vars.names[n], options="COMPRESS=LZW", overwrite=TRUE)
		
		}

		# does this need to be bias corrected here

		

		# # In this instance I have found that the Raster package has a dummed down version of the proj4 lib and it is giving strange results
		# # therefore I am using a system call from R to the GDAL library to do the reproject/resample/clipping of the data...
		# #gdalCall <- paste("gdalwarp -t_srs '", akalb, "' -tr 1000 1000 -r bilinear -srcnodata -9999 -dstnodata NoData -co 'COMPRESS=LZW' -cutline ", paste(outPath,'extent/AIEM_domain.shp', sep="")," -overwrite ", paste(outPath,"downscaled/hur_cru_ts31/10min/","hur_cru_ts31_10min_",month,"_",yearList[count],".tif",sep="")," " ,paste(outPath,"downscaled/hur_cru_ts31/1km/","hur_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), sep="")
		# # new version with out the cutline
		# gdalCall <- paste("gdalwarp -t_srs '", akalb, "' -tr 1000 1000 -r bilinear -srcnodata -9999 -dstnodata NoData -co 'COMPRESS=LZW'"," -overwrite ", paste(outPath,"downscaled/hur_cru_ts31/10min/","hur_cru_ts31_10min_",month,"_",yearList[count],".tif",sep="")," " ,paste(outPath,"downscaled/hur_cru_ts31/1km/","hur_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), sep="")
		# system(gdalCall)
		
		# # and the tas stuff
		# gdalCall <- paste("gdalwarp -t_srs '", akalb, "' -tr 1000 1000 -r bilinear -srcnodata -9999 -dstnodata NoData -co 'COMPRESS=LZW'"," -overwrite ", paste(outPath,"downscaled/tas_cru_ts31/10min/","tas_cru_ts31_10min_",month,"_",yearList[count],".tif",sep="")," " ,paste(outPath,"downscaled/tas_cru_ts31/1km/","tas_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), sep="")
		# system(gdalCall)

		# gdalCall <- paste("gdalwarp -t_srs '", akalb, "' -tr 1000 1000 -r bilinear -srcnodata NA -dstnodata NoData -co 'COMPRESS=LZW'"," -overwrite ", paste(outPath,"downscaled/vap_cru_ts31/10min/","vap_cru_ts31_10min_",month,"_",yearList[count],".tif",sep="")," " ,paste(outPath,"downscaled/vap_cru_ts31/1km/","vap_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), sep="")
		# system(gdalCall)

		# # here we are going to perform some data acrobatics to get the data to play with the existing stuff for ALF
		# # read back in the 1km hur
		# hur.mask <- raster(paste(outPath,"downscaled/hur_cru_ts31/1km/","hur_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""))
		# # crop the data to the extent of the output mask
		# hur.mask.c <- crop(hur.mask,out.mask)
		# # then set the damn extent to the out.mask
		# extent(hur.mask.c) <- extent(out.mask)
		# # now lets mask the output to match ALFRESCO stuff
		# mask(hur.mask.c, out.mask, filename=paste(outPath,"downscaled/hur_cru_ts31/1km_mask/","hur_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), options="COMPRESS=LZW", overwrite=TRUE)
		
		# # read back in the 1km vapor
		# vap.mask <- raster(paste(outPath,"downscaled/vap_cru_ts31/1km/","vap_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""))
		# # crop the data to the extent of the output mask
		# vap.mask.c <- crop(vap.mask,out.mask)
		# # then set the damn extent to the out.mask
		# extent(vap.mask.c) <- extent(out.mask)
		# # now lets mask the output to match ALFRESCO stuff
		# mask(vap.mask.c, out.mask, filename=paste(outPath,"downscaled/vap_cru_ts31/1km_mask/","vap_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), options="COMPRESS=LZW", overwrite=TRUE)

		# # read back in the 1km tas
		# tas.mask <- raster(paste(outPath,"downscaled/tas_cru_ts31/1km/","tas_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""))
		# # crop the data to the extent of the output mask
		# tas.mask.c <- crop(tas.mask,out.mask)
		# # then set the damn extent to the out.mask
		# extent(tas.mask.c) <- extent(out.mask)
		# # now lets mask the output to match ALFRESCO stuff
		# mask(tas.mask.c, out.mask, filename=paste(outPath,"downscaled/tas_cru_ts31/1km_mask/","tas_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), options="COMPRESS=LZW", overwrite=TRUE)

	}
}

d.end <- date()

print(paste("BEGIN_TIME: ",d.begin, sep=""))
print("----------------------------------")
print(paste("  END_TIME: ",d.end, sep=""))



