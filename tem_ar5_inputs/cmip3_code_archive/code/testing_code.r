# this is some testing to deal with the undoubted extent issues that will arise with the 

l <- list.files("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/cru_files", pattern="_PCLL.tif", full.names=T)

# use that list object to stack all the files to a stack
crustack <- stack(l)

# read in the extent shapefile for the interpolation analysis
extent.shape <- readShapeSpatial("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/idrisi/iem_downscale_mask_FINAL.shp")

crustack.desiredCells <- cellsFromExtent(subset(crustack,1,drop=T), extent.shape)
crustack.desiredCells.e <- crustack[crustack.desiredCells,drop=F]

ncstack.anom.spline	<- interp(x=in_xy[,1],y=in_xy[,2],z=z_in,xo=seq(min(out_xy[,1]),max(out_xy[,1]),l=ncol(cru.current)), yo=seq(min(out_xy[,2]),max(out_xy[,2]),l=nrow(cru.current),linear=F))




for(i in 1:nlayers(crustack.desiredCells.e)){
	
	#print(crustack.desiredCells.e[i])

	cru.current <- subset(crustack.desiredCells.e, i, drop=T)

	z_in <- getValues(cru.current)

	crustack.anom.spline <- interp(x=in_xy[,1],y=in_xy[,2],z=z_in,xo=seq(min(out_xy[,1]),max(out_xy[,1]),l=ncol(cru.current)), yo=seq(min(out_xy[,2]),max(out_xy[,2]),l=nrow(cru.current)),linear=F)

}

#print(crustack.desiredCells.e[i])


out_xy <- coordinates(subset(crustack.desiredCells.e,1,drop=T))

cru.current <- subset(crustack.desiredCells.e, 1, drop=T)

in_xyz <- data.frame(coordinates(subset(crustack.desiredCells.e,1,drop=T)),getValues(subset(crustack.desiredCells.e,1,drop=T)))

in_xyz <- in_xyz[in_xyz[,3] > -9999,]

crustack.anom.spline <- interp(x=in_xyz[,1], y=in_xyz[,2], z=in_xyz[,3], xo=seq(min(out_xy[,1]),max(out_xy[,1]),l=ncol(cru.current)), yo=seq(min(out_xy[,2]),max(out_xy[,2]),l=nrow(cru.current)),linear=F)

crustack.interp <- t(crustack.anom.spline$z)[,nrow(crustack.anom.spline$z):1]

nc.interp.r <- raster(crustack.interp, xmn=xmin(cru.current), xmx=xmax(cru.current), ymn=ymin(cru.current), ymx=ymax(cru.current))