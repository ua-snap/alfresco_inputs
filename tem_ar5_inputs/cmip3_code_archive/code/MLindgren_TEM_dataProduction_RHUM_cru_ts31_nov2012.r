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
l <- list.files("/big_storage/malindgren/AIEM/RHUM/hur_cru_10min_1961_1990/pcll_spline_final", pattern=".tif", full.names=T)

# use that list object to stack all the files to a stack
crustack.hur <- stack(l)

l <- list.files("/big_storage/malindgren/AIEM/RHUM/tas_cru_10min_1961_1990/pcll_spline_final", pattern=".tif", full.names=T)
crustack.tas <- stack(l)

# read in the extent shapefile for the interpolation analysis (in Pacific Centered LatLong)
extent.shape <- extent(c(181.1666, 307.3343, 41.49969, 83.33335)) # OLD >>> extent(c(166.6012,305.4721,38.023,83.4692)) #readShapeSpatial("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/WORKING_TEM/Extents/Extents_TEM_dataCreate/iem_downscale_mask_FINAL.shp")

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
print("*****  CONVERTING VAPOR PRESSURE INTO RELATIVE HUMIDITY	************")

# back to the old conversion
esa.stack = 6.112*exp(17.62*tas.stack/(243.12+tas.stack))
hur.stack <- vap.stack/esa.stack * 100


##########################################################################################################
print("***** 	FLIPPING RASTERS TO PACIFIC CENTERED LATLONG 	************")

# this little codeblock is to take the not Pacific Centered CRU 0.5 and p
# get the coordinates of the centroids and the values of the stack in a matrix
hur.xyz <- data.frame(coordinates(hur.stack), getValues(hur.stack))
tas.xyz <- data.frame(coordinates(tas.stack), getValues(tas.stack))

# now we use a bit of voo-doo to make the data sit as pacific centered lonlat
if(any(hur.xyz[,1] < 0))  hur.xyz[hur.xyz[,1]<0,1] <- hur.xyz[hur.xyz[,1]<0,1]+360
if(any(tas.xyz[,1] < 0))  tas.xyz[tas.xyz[,1]<0,1] <- tas.xyz[tas.xyz[,1]<0,1]+360

# use these as iterators to create the files you want
monList <- 1:12
yearList <- 1901:2009

# create some empty lists to grab and hold the filenaming for stacking later
hur.stack <- stack()
tas.stack <- stack()

iter=0
for(y in yearList){
	for(mon in monList){
		# this iterator will help move through the columns
		iter=iter+1

		print(iter)

		if(nchar(mon)<2){ mon=paste("0",mon,sep="")}else{mon=paste(mon,sep="")}	
		
		# this is where we make the new raster from the xyz 
		hur.flipped <- rasterFromXYZ(cbind(hur.xyz[,1:2],hur.xyz[(iter+2)]), res=c(0.5,0.5), crs=NA, digits=5)
		tas.flipped <- rasterFromXYZ(cbind(tas.xyz[,1:2],tas.xyz[(iter+2)]), res=c(0.5,0.5), crs=NA, digits=5)
		
		# here we write that converted raster object to a new tif file
		writeRaster(hur.flipped, filename=paste(outPath,"cru_historical_1901_2009/pacific/hur/","hur_total_mm_cru_TS31_pcLL_", mon,"_",y,".tif",sep=""), overwrite=T)
		writeRaster(tas.flipped, filename=paste(outPath,"cru_historical_1901_2009/pacific/tas/","tas_total_mm_cru_TS31_pcLL_", mon,"_",y,".tif",sep=""), overwrite=T)

		# keep a list of the filenames in chronological order so that we can list em using quickstack after the loop
		hur.stack <- addLayer(hur.stack, hur.flipped)
		tas.stack <- addLayer(tas.stack, tas.flipped)
		# hur.l <- append(hur.l, paste(outPath,"cru_historical_1901_2009/pacific/hur/","hur_total_mm_cru_TS31_pcLL_", mon,"_",y,".tif",sep=""), after = length(hur.l))
		# tas.l <- append(tas.l, paste(outPath,"cru_historical_1901_2009/pacific/tas/","tas_total_mm_cru_TS31_pcLL_", mon,"_",y,".tif",sep=""), after = length(tas.l))
	}
}




# this may no longer be needed <<< CHECK IT AFTER COMPLETION
# this line will re-initialize vap.stack with the newly flipped values
# hur.stack <- quickStack(hur.l)
# tas.stack <- quickStack(tas.l)

# # remove the list object(s)
# rm(hur.l) 
# rm(tas.l)

# ######################################################	END FILE LISTING LOOP  	######################################################################################################################
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ######################################################	LOOP TO FLIP TO PC Latlong 	##################################################################################################################

# count=0
# for(j in 3:ncol(hur.xyz)){
# 	print(paste("working on flipping: ", colnames(hur.xyz[j])))
# 	# here we select the XYZ for a given month/year combination in the timeseries and convert to raster
# 	hur.flipped <- rasterFromXYZ(cbind(hur.xyz[,1:2],hur.xyz[j]), res=c(0.5,0.5), crs=NA, digits=5)
# 	tas.flipped <- rasterFromXYZ(cbind(tas.xyz[,1:2],tas.xyz[j]), res=c(0.5,0.5), crs=NA, digits=5)
# 	# here we write that converted raster object to a new tif file
# 	writeRaster(hur.flipped, filename=hur.l[j-2], overwrite=T)
# 	writeRaster(tas.flipped, filename=tas.l[j-2], overwrite=T)
# }


# #####################################################  	ABOVE Not Necessary For All Runs	#########################################################################################################
# #####################################################################################################################################################################################################



# # this line will re-initialize vap.stack with the newly flipped values
# hur.stack <- quickStack(hur.l)
# tas.stack <- quickStack(tas.l)

# # remove the list object(s)
# rm(hur.l) 
# rm(tas.l)






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


print("***** 	CREATING MONTHLY CLIMATOLOGY 	************************")

# here we are going to use the subset command to grab the files for 1961-1990
#  ** these values only make sense if the 20c3m or historical timeseries dates are monthly 1901-2009
hur.stack.select <- subset(hur.stack, 721:1080, drop=T)
tas.stack.select <- subset(tas.stack, 721:1080, drop=T)

# create a new stack for the output rasters
climstack.mean.hur <- stack()
climstack.mean.tas <- stack()
# loop through the months in the years and calculate the mean monthly climatology for the period of 1961-1990
#  this climatology is created to create anomalies of the CRU TS3.1 historical data using the CRU TS3.1 1961-1990 as the baseline period
for(i in 1:12){
	monthList <- seq(i,nlayers(hur.stack.select),12) # a list of month indexes through subselected series
	if(nchar(i)<2){ month=paste("0",i,sep="")}else{month=paste(i,sep="")} # month naming convention 
	substack.hur <- subset(hur.stack.select, monthList) # subset the series
	substack.tas <- subset(tas.stack.select, monthList)

	monMean.hur <- mean(substack.hur) # run a mean to get the climatology for the 30 year period
	monMean.tas <- mean(substack.tas)

	# nameholder.hur <- paste(outPath,"hur_cru_ts31_1961_1990/","hur_cru_ts_3_10_pcLL_", month,"_1961_1990.tif",sep="")
	# nameholder.tas <- paste(outPath,"tas_cru_ts31_1961_1990/","tas_cru_ts_3_10_pcLL_", month,"_1961_1990.tif",sep="")
	
	# print(paste("writing: ",nameholder.hur,sep=""))
	writeRaster(monMean.hur, filename=paste(outPath,"hur_cru_ts31_1961_1990/","hur_cru_ts_3_10_pcLL_", month,"_1961_1990.tif",sep=""), overwrite=T)
	writeRaster(monMean.tas, filename=paste(outPath,"tas_cru_ts31_1961_1990/","tas_cru_ts_3_10_pcLL_", month,"_1961_1990.tif",sep=""), overwrite=T)

	climstack.mean.hur <- addLayer(monMean.hur)
	climstack.mean.tas <- addLayer(monMean.tas)

}

rm(nameholder.hur)
rm(nameholder.tas)

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

# # list the newly aggregated mean monthly files
l.hur <- character()
l.tas <- character()

monList <- 1:12

# I am doing this is this loop because due to the filename structure they will not "list" chronologically with list.files()

for(mon in monList){
	if(nchar(mon)<2){ mon=paste("0",mon,sep="")}else{mon=paste(mon,sep="")}	

	l.hur <- append(l.hur, paste(outPath,"hur_cru_ts31_1961_1990/","hur_cru_ts_3_10_pcLL_", mon,"_1961_1990.tif",sep=""), after = length(l.hur))
	l.tas <- append(l.tas, paste(outPath,"tas_cru_ts31_1961_1990/","tas_cru_ts_3_10_pcLL_", mon,"_1961_1990.tif",sep=""), after = length(l.tas))
}


# stack up the newly averaged files
climstack.mean.hur <- quickStack(l.hur)
climstack.mean.tas <- quickStack(l.tas)

# remove the current list of files
rm(l.hur)
rm(l.tas)

# grab the desired cells from each stack object being used for the new climatology stack
climstack.mean.hur.desiredCells <- cellsFromExtent(subset(climstack.mean.hur,1,drop=T), extent.shape)

# now extract those cells to new stacks for the new climatology stack
climstack.mean.hur.desiredCells.e <- climstack.mean.hur[climstack.mean.hur.desiredCells,drop=F]

# grab the desired cells from each stack object being used for the new climatology stack
climstack.mean.tas.desiredCells <- cellsFromExtent(subset(climstack.mean.tas,1,drop=T), extent.shape)

# now extract those cells to new stacks for the new climatology stack
climstack.mean.tas.desiredCells.e <- climstack.mean.tas[climstack.mean.tas.desiredCells,drop=F]

print("***** 	CALCULATING MONTHLY ANOMALIES 	************************")

# create anomalies 
for(i in 1:12){
	print(paste("	MONTH WORKING = ",i))
	monthList <- seq(i, nlayers(hur.stack), 12)
	if(nchar(i)<2){ month=paste("0",i,sep="")}else{month=paste(i,sep="")}
	climstack.hur.current <- subset(climstack.mean.hur.desiredCells.e, i, drop=T)
	climstack.tas.current <- subset(climstack.mean.tas.desiredCells.e, i, drop=T)
	
	count=0 # a counter to iterate through the years
	
	for(j in monthList){
		print(paste("anomalies iter ",j))
		count=count+1
		hur.stack.current <- subset(hur.stack.desiredCells.e, j, drop=T)
		tas.stack.current <- subset(tas.stack.desiredCells.e, j, drop=T)
		
		# dont know why this is necessary, but I am doing it so that the software doesnt complain about #isfactor()
		tas.stack.current@data@isfactor <- FALSE 

		anom.hur <- hur.stack.current/climstack.hur.current
		anom.tas <- tas.stack.current-climstack.tas.current

		nameholder.hur <- paste(outPath,"cru_historical_1901_2009/anomalies_pacific/hur/","hur_cru_ts_3_10_pcLL_anom_",month,"_",yearList[count],".tif",sep="")
		nameholder.tas <- paste(outPath,"cru_historical_1901_2009/anomalies_pacific/tas/","tas_cru_ts_3_10_pcLL_anom_",month,"_",yearList[count],".tif",sep="")

		writeRaster(anom.hur, filename=nameholder.hur, overwrite=T)
		writeRaster(anom.tas, filename=nameholder.tas, overwrite=T)
	}
}

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

		l.hur <- append(l.hur, paste(outPath,"cru_historical_1901_2009/anomalies_pacific/hur/","hur_cru_ts_3_10_pcLL_anom_",month,"_",i,".tif",sep=""), after = length(l.hur))
		l.tas <- append(l.tas, paste(outPath,"cru_historical_1901_2009/anomalies_pacific/tas/","tas_cru_ts_3_10_pcLL_anom_",month,"_",i,".tif",sep=""), after = length(l.tas))
	}
}

hur.stack.anom <- quickStack(l.hur)
tas.stack.anom <- quickStack(l.tas)

rm(l.hur)
rm(l.tas)

in_xy <- coordinates(subset(hur.stack.anom,1,drop=T))
out_xy <- coordinates(subset(crustack.hur.desiredCells.e,1,drop=T))


print("***** 	PERFORMING DOWNSCALING OF TIMESERIES 	******************")


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
		downscaled.month.flip.hur.r <- rasterFromXYZ(downscaled.month.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", digits=8) #downscaled.month.flip.r <- rasterFromXYZ(downscaled.month.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")


		### downscale the temperature
		hur.stack.current <- subset(tas.stack.anom, j, drop=T)

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
		downscaled.month.flip.tas.r <- rasterFromXYZ(downscaled.month.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", digits=8) #downscaled.month.flip.r <- rasterFromXYZ(downscaled.month.flip, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

		#### END TAS

		# # here we get the values of the current downscaled layer
		downscaled.month.flip.hur.r.v <- getValues(downscaled.month.flip.hur.r)
		downscaled.month.flip.tas.r.v <- getValues(downscaled.month.flip.tas.r)

		# # now we ask some questions about the values in the map.  The acceptable range is 0-100 ** although there are instances where there is >100 
		# # here we are asking it to turn all of these out of bounds values to a NoData value of -9999
		ind <- which(downscaled.month.flip.hur.r.v < 0); values(downscaled.month.flip.hur.r)[ind] <- -9999
		#ind <- which(downscaled.month.flip.tas.r.v < 0); values(downscaled.month.flip.tas.r)[ind] <- -9999

		# now we change the >100 values to 95 as per steph's okaying it
		values(downscaled.month.flip.hur.r)[which(getValues(downscaled.month.flip.hur.r) > 100)] <- 95

		print("***** 	COVNERTING RELATIVE HUMIDITY BACK TO VAPOR PRESSURE 	******************")

		# convert back to vapor pressure
		esa = 6.112*exp(17.62*downscaled.month.flip.tas.r/(243.12+downscaled.month.flip.tas.r))
		vapor.ts31.10min = (downscaled.month.flip.hur.r*esa)/100


		# here we write out the greenwich-centered version of the map
		writeRaster(trim(downscaled.month.flip.hur.r), filename=paste(outPath,"downscaled/hur_cru_ts31/10min/","hur_cru_ts31_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)
		writeRaster(trim(downscaled.month.flip.tas.r), filename=paste(outPath,"downscaled/tas_cru_ts31/10min/","tas_cru_ts31_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)
		writeRaster(trim(vapor.ts31.10min), filename=paste(outPath,"downscaled/vap_cru_ts31/10min/","vap_cru_ts31_10min_",month,"_",yearList[count],".tif",sep=""), overwrite=T)


		# this is the proj4 string to reproject into Alaska Albers 
		newProj <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +nadgrids=@alaska +datum=NAD83 +units=m +no_defs"

		# In this instance I have found that the Raster package has a dummed down version of the proj4 lib and it is giving strange results
		# therefore I am using a system call from R to the GDAL library to do the reproject/resample/clipping of the data...
		#gdalCall <- paste("gdalwarp -t_srs '", newProj, "' -tr 1000 1000 -r bilinear -srcnodata -9999 -dstnodata NoData -co 'COMPRESS=LZW' -cutline ", paste(outPath,'extent/AIEM_domain.shp', sep="")," -overwrite ", paste(outPath,"downscaled/hur_cru_ts31/10min/","hur_cru_ts31_10min_",month,"_",yearList[count],".tif",sep="")," " ,paste(outPath,"downscaled/hur_cru_ts31/1km/","hur_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), sep="")
		# new version with out the cutline
		paste(outPath,"downscaled/hur_cru_ts31/10min/","hur_cru_ts31_10min_",month,"_",yearList[count],".tif",sep="")

		gdalCall <- paste("gdalwarp -t_srs '", newProj, "' -tr 1000 1000 -r bilinear -srcnodata -9999 -dstnodata NoData -co 'COMPRESS=LZW'"," -overwrite ", paste(outPath,"downscaled/hur_cru_ts31/10min/","hur_cru_ts31_10min_",month,"_",yearList[count],".tif",sep="")," " ,paste(outPath,"downscaled/hur_cru_ts31/1km/","hur_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), sep="")

		# this is how we issue a system call in R.
		system(gdalCall)
		
		# and the tas stuff
		gdalCall <- paste("gdalwarp -t_srs '", newProj, "' -tr 1000 1000 -r bilinear -srcnodata -9999 -dstnodata NoData -co 'COMPRESS=LZW'"," -overwrite ", paste(outPath,"downscaled/tas_cru_ts31/10min/","tas_cru_ts31_10min_",month,"_",yearList[count],".tif",sep="")," " ,paste(outPath,"downscaled/tas_cru_ts31/1km/","tas_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), sep="")

		system(gdalCall)

		gdalCall <- paste("gdalwarp -t_srs '", newProj, "' -tr 1000 1000 -r bilinear -srcnodata NA -dstnodata NoData -co 'COMPRESS=LZW'"," -overwrite ", paste(outPath,"downscaled/vap_cru_ts31/10min/","vap_cru_ts31_10min_",month,"_",yearList[count],".tif",sep="")," " ,paste(outPath,"downscaled/vap_cru_ts31/1km/","vap_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), sep="")

		system(gdalCall)

		# read back in the 1km hur
		hur.mask <- raster(paste(outPath,"downscaled/hur_cru_ts31/1km/","hur_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""))
		# read back in the 1km vapor
		vap.mask <- raster(paste(outPath,"downscaled/vap_cru_ts31/1km/","vap_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""))
		# read back in the 1km tas
		tas.mask <- raster(paste(outPath,"downscaled/tas_cru_ts31/1km/","tas_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""))

		# and now we mask them and write 'em out
		cropped <- crop(hur.mask,out.mask)
		mask(hur.mask, out.mask, filename=paste(outPath,"downscaled/hur_cru_ts31/1km_mask/","hur_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), options="COMPRESS=LZW")
		mask(vap.mask, out.mask, filename=paste(outPath,"downscaled/vap_cru_ts31/1km_mask/","vap_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), options="COMPRESS=LZW")
		mask(tas.mask, out.mask, filename=paste(outPath,"downscaled/tas_cru_ts31/1km_mask/","tas_cru_ts31_1km_",month,"_",yearList[count],".tif",sep=""), options="COMPRESS=LZW")
	}
}

d.end <- date()

print(paste("BEGIN_TIME: ",d.begin, sep=""))
print("----------------------------------")
print(paste("  END_TIME: ",d.end, sep=""))
