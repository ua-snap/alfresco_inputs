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
setOptions(tmpdir="/tmp/")
setOptions(progress="text")
setOptions(timer=TRUE)

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


#list the cru data 

l <- list.files("/big_storage/malindgren/AIEM/RHUM/hur_cru_10min_1961_1990", pattern="_PCLL.tif", full.names=T)
crustack <- quickStack(l)
# read in the extent shapefile for the interpolation analysis (in Pacific Centered LatLong)
extent.shape <- extent(c(166.6012,305.4721,38.023,83.4692))

crustack.desiredCells <- cellsFromExtent(subset(crustack,1,drop=T), extent.shape)
crustack.desiredCells.e <- crustack[crustack.desiredCells,drop=F]


out_xy <- coordinates(crustack.desiredCells.e)

for(j in 1:nlayers(crustack.desiredCells.e)){
	print(j)
	cur <- subset(crustack.desiredCells.e, j, drop=T)

	in_xyz <- cbind(out_xy, getValues(cur))
	in_xyz = in_xyz[is.na(in_xyz[,3]) == FALSE,]

	#in_xyz <- na.omit(data.frame(in_xy, getValues(hur.stack.current)))

	hur.stack.anom.spline <- interp(x=in_xyz[,1],y=in_xyz[,2],z=in_xyz[,3],xo=seq(min(out_xy[,1]),max(out_xy[,1]),l=ncol(cur)), yo=seq(min(out_xy[,2]),max(out_xy[,2]),l=nrow(cur)))

	# transpose the data here 
	# check this result!
	nc.interp <- t(hur.stack.anom.spline$z)[,nrow(hur.stack.anom.spline$z):1]

	# rasterize it
	nc.interp.r <- raster(nc.interp, xmn=xmin(cru.current), xmx=xmax(cru.current), ymn=ymin(cru.current), ymx=ymax(cru.current))
	writeRaster(nc.interp, filename=paste("/big_storage/malindgren/AIEM/RHUM/hur_cru_10min_1961_1990/cru_10min_splined_extent/",crustack@layernames[j],".tif",sep=""))
}
