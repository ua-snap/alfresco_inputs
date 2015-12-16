# bring in the needed libraries
require(aspace)
require(raster)
require(sp)
require(maptools)
require(rgdal)


rsds <- brick("/Data/Base_Data/Climate/World/GCM_raw/rsds/20c3m/pcmdi.ipcc4.mpi_echam5.20c3m.run1.monthly.rsds_A1.nc")

# change the extent of this layer into something that actually makes sense in the pacific centered latlong system
extent(rsds) <- c(0,360,-89.50434,89.50434)

rsds.rot <- rotate(rsds)

rsds.rot.select <- subset(rsds.rot, 1:12, drop=TRUE)

latlong <- coordinates(rsds.rot)

dateList <- rsds.rot.select@layernames
ordList <- numeric()

for(i in 1:length(dateList)){
	year <- substr(dateList[i],2,5)
	mon <- substr(dateList[i],7,8)
	day <- substr(dateList[i],10,11)

	newDate <- paste(mon,day,year, sep="/")
	#x <- newDate
	# convert to day of year (Julian date) -- use POSIXlt
	ordCur <- strptime(newDate, "%m/%d/%Y")$yday+1

	#if(nchar(ordCur)==1){ ordCur <- paste("00",ordCur,sep="") } else if(nchar(ordCur)==2){ ordCur <- paste("0",ordCur,sep="") } else { ordCur <- ordCur }

	print(paste("working on Ordinal Day: ", ordCur, sep=""))

	#ordList <- append(ordList, ordCur, after=length(ordList))

	lat2rad <- as_radians(latlong[,2])

	out.matrix <- calcRa(ordCur,lat2rad,nrow(rsds.rot),ncol(rsds.rot))

	out.matrix.t <- t(out.matrix)

	out.rast <- raster(out.matrix, xmn=xmin(rsds.rot), xmx=xmax(rsds.rot), ymn=ymin(rsds.rot), ymx=ymax(rsds.rot))

	# writeout the new raster file 
	writeRaster(out.rast, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/rad_calculation/","rad_ar4_mpi_echam5_20c3m_",mon,"_",year,"_NEW.tif", sep=""), overwrite=TRUE)
}



test <- t(out.matrix)[,nrow(out.matrix):1]

test <- out.matrix[,nrow(out.matrix):1]

in_r <-raster("/Data/Base_Data/Climate/AK_CAN_2km/AK_CAN_2km/cru_TS31/historical/tas/tas_mean_C_cru_TS31_12_2009.tif")

alb<-coordinates(in_r)
alb_y <-alb[,2]
alb_y_rad<-as_radians(alb_y)
new <- calcRa(227,alb_y_rad,nrow(in_r),ncol(in_r))
image(new)

