library(raster)
library(akima)
library(sp)
library(rgdal)

# read in the netCDF to a stack
ncstack <- stack("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.mpi_echam5.sresa1b.run1.monthly.hur_A1_2001-2050.nc")

# list the tiff files from the directory containing the CRU climatologies for the variable being downscaled
l <- list.files("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/extracted_files", pattern=".tif$", full.names=T)

# use that list object to stack all the files to a stack
crustack <- stack(l)

# now lets rotate the cru data to fit the pacific centered latlong
if(any(xy[,1]>180)) xy[,1][xy[,1]>180] <- xy[,1][xy[,1]>180]-360
# cru stuff
if(any(cru.out.xy[,1]>180)) cru.out.xy[,1][cru.out.xy[,1]>180] <- cru.out.xy[,1][cru.out.xy[,1]>180]-360


cru <- raster(crustack, layer=1)
cru.out.xy <- coordinates(cru)

# now that we have something with an accurate ref system and we know that the ref system of the ncstack is WGS1984 Greenwich
# lets give the data a projection
#projection(ncstack) <- projection(cru)

#######

for(i in 1:12){ # this is a loop that creates an index of the monthly files I want to grab from the ncstack on monthly basis
	monthList <- seq(i, nlayers(ncstack), 12)
		
	cru <- raster(crustack, layer=i)
	cru.out.xy <- coordinates(cru)

	for(j in monthList){ # this loop iterates through the index of all of the * month in the series
		print(j)
		
		# grab the jth raster from the stack
		r <- raster(ncstack, layer=j)
		
		# rotate the raster to the greenwich centering Latlong
		r.rot <- rotate(r)

		# Change the projection system
		projection(r.rot) <- projection(cru)

		#crop that little bugger to the extent of the CRU raster
		r.crop <- crop(r.rot, extent(cru))

		xy <- cbind(coordinates(r.crop), getValues(r.crop))
		colnames(xy) <- c("lon","lat","rhum")

		# flip the outputs to the standard greenwich latlong like the CRU data
		#if(any(xy[,1]>180)) xy[,1][xy[,1]>180] <- xy[,1][xy[,1]>180]-360# if(any(xy[,1]>180)) xy[,1][xy[,1]>180] <- xy[,1][xy[,1]>180]-360

		# akima interpolate that
		r.i <- interp(xy[,1], xy[,2], xy[,3], xo=seq(min(cru.out.xy[,1]), max(cru.out.xy[,1]), l=nrow(cru)), yo=seq(min(cru.out.xy[,2]),max(cru.out.xy[,2]), l=ncol(cru)), linear=TRUE)

		#flip it!  This comes from matt leonowicz
		flipped <- t(r.i$z)[,nrow(r.i$z):1]

		# make a new raster layer object with the bounding coordinates from the cru 10min
		out <- raster(flipped, xmn=xmin(cru), xmx=xmax(cru), ymn=ymin(cru), ymx=ymax(cru))

		# Give the proper projection to flipped raster
		projection(out) <- projection(cru)

		# write the little bugger to a file
		writeRaster(out, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/downscaled/","hur_cru_10min_alf_cccma_",j,".tif", sep=""), overwrite=T) #pr_total_mm_alf_cru_TS31_01_1901
	}
}


# ######

# for(i in 1:nlayers(ncstack)){
# 	count=count+1

# 	# grab the first layer in the timeseries through this selection procedure
# 	r <- raster(ncstack, layer=i)

# 	# here we get the "xy" from each cell in the input low res raster
# 	#xy <- data.frame(xyFromCell(r, 1:ncell(r)))
# 	xy <- data.frame(coordinates(r), getValues(r))
# 	colnames(xy) <- c("x","y","rhum")

# 	# here we ask for the values of that low res raster
# 	#r.v <- getValues(r)

# 	#-----------------------------------------------------------------------------
# 	# here I am creating a dummy map to interpolate to for testing.
# 	# disagg
# 	#cru <- disaggregate(r, 10, method='bilinear')

# 	# get the xy coords of the cells
# 	#cru.out.xy <- data.frame(xyFromCell(cru, 1:ncell(cru)))

# 	#-----------------------------------------------------------------------------
# 	# akima interpolate that
# 	r.i <- interp(xy[,1], xy[,2], xy[,3], xo=seq(min(cru.out.xy[,1]), max(cru.out.xy[,1]), l=nrow(cru)), yo=seq(min(cru.out.xy[,2]),max(cru.out.xy[,2]), l=ncol(cru)), linear=TRUE)

# 	# now we turn the little bugger into a raster file that can be brought into arcgis or another GIS
# 	new <- raster(r.i$z)

# 	writeRaster(new, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/downscaled/","rhum_10min_cccma_a1b_",count,".tif",sep=""), overwrite=T)
# }

# #this is the code from the actual docs for testing

# interp(x, y, z, xo=seq(min(x), max(x), length = 40), yo=seq(min(y), max(y), length = 40), linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL, ncp = NULL)


# tas.i <- interp(xy[,1], xy[,2], tas.v, xo=seq(min(xy.hi[,1]), max(xy.hi[,1]), l=ncol(tas.hi)), yo=seq(min(xy.hi[,2]),max(xy.hi[,2]), l=nrow(tas.hi)), linear=TRUE)