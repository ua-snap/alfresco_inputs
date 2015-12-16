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
setOptions(tmpdir="/tmp/", progress="text", timer=TRUE)

# here we bring in some functions to be used in the processing...
source("/big_storage/malindgren/AIEM/functions/CalcAnom.r")
source("/big_storage/malindgren/AIEM/functions/ListFiles.r")

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
# vap.stack <- brick("/Data/Base_Data/Climate/World/CRU_grids/CRU_TS31/cru_ts_3_10.1901.2009.vap.dat.nc")
# tas.stack <- brick("/Data/Base_Data/Climate/World/CRU_grids/CRU_TS31/cru_ts_3_10.1901.2009.tmp.dat.nc")

# list the tiff files from the directory containing the CRU climatologies for the variable being downscaled
#  This line assumes some things about the input data, as in I have already converted these GLOBAL data into
#  Pacific centered Lat long. 
# l <- list.files("/big_storage/malindgren/AIEM/RHUM/hur_cru_10min_1961_1990/pcll_spline_final", pattern=".tif", full.names=T)

# use that list object to stack all the files to a stack
crustack.hur <- stack(list.files("/big_storage/malindgren/AIEM/RHUM/climatologies/hur/hur_cru_10min_1961_1990/pcll_spline_final", pattern=".tif", full.names=T))

# l <- list.files("/big_storage/malindgren/AIEM/RHUM/tas_cru_10min_1961_1990/pcll_spline_final", pattern=".tif", full.names=T)
crustack.tas <- stack(list.files("/big_storage/malindgren/AIEM/RHUM/climatologies/tas/tas_cru_10min_1961_1990/pcll_spline_final", pattern=".tif", full.names=T))

# read in the extent shapefile for the interpolation analysis (in Pacific Centered LatLong)
extent.shape <- extent(c(160, 310, 30, 84)) # OLD >>> extent(c(166.6012,305.4721,38.023,83.4692)) #readShapeSpatial("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/WORKING_TEM/Extents/Extents_TEM_dataCreate/iem_downscale_mask_FINAL.shp")

# this is a list of the years that is used in creating the output naming convention
yearList <- 1901:2009 

# cru 10 min akalb clipped template
template10min <- raster(raster("/big_storage/malindgren/AIEM/RHUM/templates/template_raster_10min_akalb.tif"))
template10min.wgs <- raster(raster("/big_storage/malindgren/AIEM/RHUM/templates/template_raster_10min_wgs84.tif"))


# list the uninterpolated anomalies (PCLL)
hur.l <- listem("/big_storage/malindgren/AIEM/RHUM/anomalies/hur/ts31/","hur_cru_ts_3_10_pcLL_anom_",1901,2009,".tif")
tas.l <- listem("/big_storage/malindgren/AIEM/RHUM/anomalies/tas/ts31/","tas_cru_ts_3_10_pcLL_anom_",1901,2009,".tif")
# stack em
hur.stack <- stack(hur.l)
tas.stack <- stack(tas.l)
# remove the lists
rm(hur.l); rm(tas.l)

# read in the CRU TS20 climatology data (PCLL)
crustack.hur <- stack(list.files("/big_storage/malindgren/AIEM/RHUM/climatologies/hur/hur_cru_10min_1961_1990/pcll_spline_final", pattern=".tif", full.names=T))
crustack.tas <- stack(list.files("/big_storage/malindgren/AIEM/RHUM/climatologies/tas/tas_cru_10min_1961_1990/pcll_spline_final", pattern=".tif", full.names=T))


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
# print("*****  CONVERTING VAPOR PRESSURE INTO RELATIVE HUMIDITY	************")

# back to the old conversion
# esa.stack = 6.112*exp(17.62*tas.stack/(243.12+tas.stack))
# hur.stack <- vap.stack/esa.stack * 100

# # get the values in a matrix of all of the stacked data
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


# # we are using this to read in the data
# iter=0
# for(y in yearList){
# 	for(mon in monList){
# 		# this iterator will help move through the columns
# 		iter=iter+1

# 		print(iter)

# 		if(nchar(mon)<2){ mon=paste("0",mon,sep="")}else{mon=paste(mon,sep="")}	
		
# 		# # this is where we make the new raster from the xyz 
# 		# hur.flipped <- rasterFromXYZ(cbind(hur.xyz[,1:2],hur.xyz[(iter+2)]), res=c(0.5,0.5), crs=NA, digits=5)
# 		# tas.flipped <- rasterFromXYZ(cbind(tas.xyz[,1:2],tas.xyz[(iter+2)]), res=c(0.5,0.5), crs=NA, digits=5)
		
# 		# # here we write that converted raster object to a new tif file
# 		# writeRaster(hur.flipped, filename=paste(outPath,"cru_historical_1901_2009/pacific/hur/","hur_total_mm_cru_TS31_pcLL_", mon,"_",y,".tif",sep=""), overwrite=T)
# 		# writeRaster(tas.flipped, filename=paste(outPath,"cru_historical_1901_2009/pacific/tas/","tas_total_mm_cru_TS31_pcLL_", mon,"_",y,".tif",sep=""), overwrite=T)

# 		# keep a list of the filenames in chronological order so that we can list em using quickstack after the loop
# 		# hur.stack <- addLayer(hur.stack, hur.flipped)
# 		# tas.stack <- addLayer(tas.stack, tas.flipped)
# 		hur.l <- append(hur.l, paste(outPath,"cru_historical_1901_2009/pacific/hur/","hur_total_mm_cru_TS31_pcLL_", mon,"_",y,".tif",sep=""), after = length(hur.l))
# 		tas.l <- append(tas.l, paste(outPath,"cru_historical_1901_2009/pacific/tas/","tas_total_mm_cru_TS31_pcLL_", mon,"_",y,".tif",sep=""), after = length(tas.l))
# 	}
# }


print("***** 	SELECTING RASTER CELLS USED IN ANALYSIS 	************")
# Here we extract the cells of interest for this calculation from the 2 bricks 
hur.stack.desiredCells <- cellsFromExtent(subset(hur.stack,1,drop=T), extent.shape)
tas.stack.desiredCells <- cellsFromExtent(subset(tas.stack,1,drop=T), extent.shape)
crustack.hur.desiredCells <- cellsFromExtent(subset(crustack.hur,1,drop=T), extent.shape)
crustack.tas.desiredCells <- cellsFromExtent(subset(crustack.tas,1,drop=T), extent.shape)

hur.stack.desiredCells.e <- hur.stack[hur.stack.desiredCells,drop=F]
tas.stack.desiredCells.e <- tas.stack[tas.stack.desiredCells,drop=F]
crustack.hur.desiredCells.e <- crustack.hur[crustack.hur.desiredCells,drop=F]
crustack.tas.desiredCells.e <- crustack.tas[crustack.tas.desiredCells,drop=F]


# # I think that I should re-project the ts20 data
# # write out the spatially subset data here:
# for(i in 1:nlayers(crustack.hur.desiredCells.e)){
# 	hur.sub <- subset(crustack.hur.desiredCells.e,i,T)
# 	hur.sub.outname <- paste("/big_storage/malindgren/AIEM/RHUM/climatologies/hur/hur_cru_10min_1961_1990/pcll_windowed_temporary/", hur.sub@layernames, "_window.tif",sep="")
# 	hur.reproj.outname <- paste("/big_storage/malindgren/AIEM/RHUM/climatologies/hur/hur_cru_10min_1961_1990/akalbers_windowed/", hur.sub@layernames, "_akalb.tif",sep="")
# 	writeRaster(rotate(hur.sub), filename=hur.sub.outname, options="COMPRESS=LZW")
# 	gdalCall <- paste("gdalwarp -t_srs '", newProj, "' -tr 20000 20000 -r bilinear -srcnodata -9999 -dstnodata NoData -co 'COMPRESS=LZW'"," -overwrite ",hur.sub.outname," " ,hur.reproj.outname, sep="")
# 	system(gdalCall)

# 	tas.sub <- subset(crustack.tas.desiredCells.e,i,T)
# 	tas.sub.outname <- paste("/big_storage/malindgren/AIEM/RHUM/climatologies/tas/tas_cru_10min_1961_1990/pcll_windowed_temporary/", tas.sub@layernames, "_window.tif",sep="")
# 	tas.reproj.outname <- paste("/big_storage/malindgren/AIEM/RHUM/climatologies/tas/tas_cru_10min_1961_1990/akalbers_windowed/", tas.sub@layernames, "_akalb.tif",sep="")
# 	writeRaster(rotate(tas.sub), filename=, options="COMPRESS=LZW")
# 	gdalCall <- paste("gdalwarp -t_srs '", newProj, "' -tr 20000 20000 -r bilinear -srcnodata -9999 -dstnodata NoData -co 'COMPRESS=LZW'"," -overwrite ", tas.sub.outname," " ,tas.reproj.outname, sep="")
# 	system(gdalCall)
# }

# crustack.tas.desiredCells.e <- 

# print("***** 	CREATING MONTHLY CLIMATOLOGY & CALCULATING ANOMALIES 	************************")

# # here we use the calcAnom() function I wrote and sourced in above
# hur.stack.anom <- calcAnom(hur.stack.desiredCells.e,721,1080,absolute=FALSE)
# tas.stack.anom <- calcAnom(tas.stack.desiredCells.e,721,1080,absolute=TRUE)

# # some grids we will use in the interpolation step
in_xy <- coordinates(subset(hur.stack.desiredCells.e,1,drop=T))
out_xy <- coordinates(subset(crustack.hur.desiredCells.e,1,drop=T))
out_xy_1km <- coordinates(out.mask)

print("***** 	PERFORMING DOWNSCALING OF TIMESERIES 	******************")

for(i in 1:12){ # this is a loop that creates an index of the monthly files I want to grab from the vap.stack on monthly basis
	d.begin <- proc.time()
	monthList <- seq(i, nlayers(hur.stack.desiredCells.e), 12)
	
	if(nchar(i)<2){ month=paste("0",i,sep="")}else{month=paste(i,sep="")}

	print(month)

	cru.hur.current <- subset(crustack.hur.desiredCells.e, i, drop=T)
	cru.tas.current <- subset(crustack.tas.desiredCells.e, i, drop=T)

	count=0 # a little counter

	for(j in monthList){
		count=count+1
	
		# downscale the relative humidity	
		hur.stack.current <- subset(hur.stack.desiredCells.e, j, drop=T)

		in_xyz <- na.omit(data.frame(in_xy, getValues(hur.stack.current)))

		hur.stack.anom.spline <- interp(x=in_xyz[,1],y=in_xyz[,2],z=in_xyz[,3],xo=seq(min(out_xy[,1]),max(out_xy[,1]),l=ncol(cru.hur.current)), yo=seq(min(out_xy[,2]),max(out_xy[,2]),l=nrow(cru.hur.current)))

		# transpose the data here 
		# check this result!
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
		downscaled.month.flip.hur.r <- rasterFromXYZ(downscaled.month.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", digits=0, res=res(template10min.wgs)) #downscaled.month.flip.r <- rasterFromXYZ(downscaled.month.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
		
		# pts.hur.spdf <- SpatialPointsDataFrame(SpatialPoints(downscaled.month.flip[,1:2], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")), data=as.data.frame(downscaled.month.flip))
		# downscaled.month.flip.hur.r <- rasterize(pts.hur.spdf,template10min.wgs, field="V3")

		### downscale the temperature
		tas.stack.current <- subset(tas.stack.desiredCells.e, j, drop=T)

		in_xyz <- na.omit(data.frame(in_xy, getValues(tas.stack.current)))

		tas.stack.anom.spline <- interp(x=in_xyz[,1],y=in_xyz[,2],z=in_xyz[,3],xo=seq(min(out_xy[,1]),max(out_xy[,1]),l=ncol(cru.tas.current)), yo=seq(min(out_xy[,2]),max(out_xy[,2]),l=nrow(cru.tas.current)))

		# transpose the data here 
		# check this result!
		nc.interp <- matrix(tas.stack.anom.spline$z,nrow=nrow(cru.tas.current),ncol=ncol(cru.tas.current),byrow=T)[nrow(cru.tas.current):1,]

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
		downscaled.month.flip.tas.r <- rasterFromXYZ(downscaled.month.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", digits=0, res=res(template10min.wgs)) #downscaled.month.flip.r <- rasterFromXYZ(downscaled.month.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
		
		# make the shapefile
		#pts.tas <- SpatialPoints(downscaled.month.flip[,1:2], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
		# pts.tas.spdf <- SpatialPointsDataFrame(SpatialPoints(downscaled.month.flip[,1:2], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")), data=as.data.frame(downscaled.month.flip))
		# downscaled.month.flip.tas.r <- rasterize(pts.tas.spdf,template10min.wgs, field="V3", update=T)
		
		#### END TAS

		# # here we get the values of the current downscaled layer
		downscaled.month.flip.hur.r.v <- getValues(downscaled.month.flip.hur.r)
		#downscaled.month.flip.tas.r.v <- getValues(downscaled.month.flip.tas.r)

		# # now we ask some questions about the values in the map.  The acceptable range is 0-100 ** although there are instances where there is >100 
		# # here we are asking it to turn all of these out of bounds values to a NoData value of -9999
		# ind <- which(downscaled.month.flip.hur.r.v < 0); values(downscaled.month.flip.hur.r)[ind] <- -9999
		# ind <- which(downscaled.month.flip.tas.r.v < 0); values(downscaled.month.flip.tas.r)[ind] <- -9999

		# now we change the >100 values to 95 as per steph's okaying it
		values(downscaled.month.flip.hur.r)[which(values(downscaled.month.flip.hur.r) > 100)] <- 95


		print("***** 	COVNERTING RELATIVE HUMIDITY BACK TO VAPOR PRESSURE 	******************")

		# make the rasters into matrices
		tas.mat <- getValues(downscaled.month.flip.tas.r) 
		hur.mat <- getValues(downscaled.month.flip.hur.r)

		# new rasters to hold the output based on the information from the current set
		vapor.ts31.10min <- raster(downscaled.month.flip.tas.r)

		# convert back to vapor pressure
		esa = 6.112*exp(17.62*tas.mat/(243.12+tas.mat))
		vapor.ts31.10min.hold = (hur.mat*esa)/100
		# set the converted values into a vapor pressure variable
		vapor.ts31.10min <- setValues(vapor.ts31.10min,vapor.ts31.10min.hold)

		# # HERE WE DO ANOTHER INTERPOLATION TO THE 1KM RESOLUTION
		# # do the hur interp here
		# in_xyz <- na.omit(data.frame(coordinates(downscaled.month.akalb.hur), getValues(downscaled.month.akalb.hur)))
		# hur.stack.anom.spline <- interp(x=in_xyz[,1],y=in_xyz[,2],z=in_xyz[,3],xo=seq(min(out_xy[,1]),max(out_xy[,1]),l=ncol(out.mask)), yo=seq(min(out_xy_1km[,2]),max(out_xy_1km[,2]),l=nrow(out.mask)))
		
		# #hur.stack.anom.spline <- interp(x=in_xyz[,1],y=in_xyz[,2],z=in_xyz[,3],xo=seq(min(out_xy[,1]),max(out_xy[,1]),l=ncol(cru.hur.current)), yo=seq(min(out_xy[,2]),max(out_xy[,2]),l=nrow(cru.hur.current)))

		# nc.interp <- matrix(hur.stack.anom.spline$z,nrow=nrow(out.mask),ncol=ncol(out.mask),byrow=T)[nrow(out.mask):1,]
		# nc.interp.r <- raster(nc.interp, xmn=xmin(out.mask), xmx=xmax(out.mask), ymn=ymin(out.mask), ymx=ymax(out.mask))
		
		# nc.interp.r <- raster(matrix(hur.stack.anom.spline$z, nrow=nrow(out.mask), ncol=ncol(out.mask)), xmn=xmin(out.mask), xmx=xmax(out.mask), ymn=ymin(out.mask), ymx=ymax(out.mask))


		# do the vap interp here


		# here we write out the greenwich-centered version of the map
		writeRaster(trim(downscaled.month.flip.hur.r), filename=paste(outPath,"downscaled/hur_cru_ts31/10min/","hur_cru_ts31_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)
		writeRaster(trim(downscaled.month.flip.tas.r), filename=paste(outPath,"downscaled/tas_cru_ts31/10min/","tas_cru_ts31_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)
		writeRaster(trim(vapor.ts31.10min), filename=paste(outPath,"downscaled/vap_cru_ts31/10min/","vap_cru_ts31_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)

		# this is the proj4 string to reproject into Alaska Albers 
		newProj <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +nadgrids=@alaska +datum=NAD83 +units=m +no_defs"

		# do the reproject here
		print("reproj 2 akalb")
		downscaled.month.akalb.hur <- projectRaster(trim(downscaled.month.flip.hur.r), res=20000, crs=newProj, method="ngb", filename=paste(outPath,"downscaled/hur_cru_ts31/10min/","hur_cru_ts31_10min_akalb_",month,"_",yearList[count],".tif",sep=""), options="COMPRESS=LZW", overwrite=T)
		#downscaled.month.akalb.tas <- projectRaster(trim(downscaled.month.flip.tas.r), res=20000, crs=newProj, method="ngb", filename=paste(outPath,"downscaled/tas_cru_ts31/10min/","hur_cru_ts31_10min_akalb_",month,"_",yearList[count],".tif",sep=""), options="COMPRESS=LZW", overwrite=T)
		vapor.ts31.10min.akalb <- projectRaster(trim(vapor.ts31.10min), res=20000, crs=newProj, method="ngb", filename=paste(outPath,"downscaled/vap_cru_ts31/10min/","vap_cru_ts31_10min_akalb_",month,"_",yearList[count],".tif",sep=""), options="COMPRESS=LZW", overwrite=T)

		# now we need to resample the data to the 1km resolution
		print("resample to 1km")
		downscaled.month.akalb.hur.1km <- resample(downscaled.month.akalb.hur,out.mask,method="ngb") #, filename=paste(outPath,"downscaled/hur_cru_ts31/1km/","hur_cru_ts31_10min_akalb_",month,"_",yearList[count],sep=""), overwrite=T, options="COMPRESS=LZW"
		# downscaled.month.akalb.tas.1km <- resample(downscaled.month.akalb.tas,out.mask,method="ngb")
		vapor.ts31.1km.akalb <- resample(vapor.ts31.10min.akalb,out.mask,method="ngb")

		# and now we mask them and write 'em out
		print("masking")
		#cropped <- crop(hur.mask,out.mask)
		mask(downscaled.month.akalb.hur.1km, out.mask, filename=paste(outPath,"downscaled/hur_cru_ts31/1km_mask/","hur_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), overwrite=T, options="COMPRESS=LZW")
		# mask(downscaled.month.akalb.tas.1km, out.mask, filename=paste(outPath,"downscaled/tas_cru_ts31/1km_mask/","tas_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), overwrite=T, options="COMPRESS=LZW")
		mask(vapor.ts31.1km.akalb, out.mask, filename=paste(outPath,"downscaled/vap_cru_ts31/1km_mask/","vap_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), overwrite=T, options="COMPRESS=LZW")
		d.end <- date()
		print(paste("calc time: ",d.begin - proc.time(),sep=""))
		
	}
}

d.end <- date()

print(paste("BEGIN_TIME: ",d.begin, sep=""))
print("----------------------------------")
print(paste("  END_TIME: ",d.end, sep=""))
